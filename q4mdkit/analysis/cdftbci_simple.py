from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import mdtraj as md
import numpy as np
import yaml

from .cdftb import (
    BOHR_PER_ANG,
    CDFTBConfig,
    FragmentConfig,
    extract_atom_types_from_hsd,
    iter_qm_coordinates,
    load_mm_indices,
    load_pccharges,
    load_qm_indices,
    parse_atom_range,
    run_dftb_in_subprocess,
)
from .cdftbci import (
    DEFAULT_RETRY_ATTEMPTS,
    UnrestrictedOrbitalData,
    compute_cdftbci_for_frame,
    compute_transfer_integral_unrestricted,
    run_cdftb_with_retry,
    setup_work_directory,
)
from .odin_overlap import compute_cross_overlap_odin


@dataclass
class SimpleCDFTBCIConfig(CDFTBConfig):
    ci_N_A: float = 101.0
    ci_N_B: float = 101.0
    ci_M_A: float = 0.0
    ci_M_B: float = 0.0
    work_directory: str = "work"
    output_file: Optional[Path] = None
    use_previous_charges: bool = False
    retry_attempts: Optional[List[Dict[str, Any]]] = None
    odin_executable: Optional[str] = None
    odin_sk_prefix: str = ""
    odin_sk_separator: str = "-"
    odin_sk_suffix: str = ".skf"
    odin_lmax: Dict[str, int] = field(default_factory=dict)
    odin_work_subdir: str = "odin_work"
    odin_keep_files: bool = False


@dataclass
class DeterminantResult:
    alpha: float
    beta: float
    total: float
    sign: int


@dataclass
class GaugeState:
    occ_alpha: np.ndarray
    occ_beta: np.ndarray


def _occ_block(coefficients: np.ndarray, n_occ: int) -> np.ndarray:
    if n_occ <= 0:
        return np.zeros((coefficients.shape[0], 0), dtype=coefficients.dtype)
    return coefficients[:, :n_occ]


def _apply_state_sign(orbital_data: UnrestrictedOrbitalData, sign: int) -> GaugeState:
    occ_alpha = _occ_block(orbital_data.C_alpha, orbital_data.n_alpha).copy()
    occ_beta = _occ_block(orbital_data.C_beta, orbital_data.n_beta).copy()

    if sign < 0:
        if occ_alpha.shape[1] > 0:
            occ_alpha[:, 0] *= -1.0
        elif occ_beta.shape[1] > 0:
            occ_beta[:, 0] *= -1.0

    return GaugeState(occ_alpha=occ_alpha, occ_beta=occ_beta)


def _determinant(matrix: np.ndarray) -> float:
    if matrix.shape[0] == 0:
        return 1.0
    return float(np.linalg.det(matrix))


def _compute_state_determinant(
    previous: GaugeState,
    current: UnrestrictedOrbitalData,
    s_ao_cross: np.ndarray,
) -> DeterminantResult:
    current_alpha = _occ_block(current.C_alpha, current.n_alpha)
    current_beta = _occ_block(current.C_beta, current.n_beta)

    if previous.occ_alpha.shape[1] != current_alpha.shape[1]:
        raise ValueError("Alpha occupation count changed between frames")
    if previous.occ_beta.shape[1] != current_beta.shape[1]:
        raise ValueError("Beta occupation count changed between frames")

    overlap_alpha = previous.occ_alpha.T @ s_ao_cross @ current_alpha
    overlap_beta = previous.occ_beta.T @ s_ao_cross @ current_beta

    det_alpha = _determinant(overlap_alpha)
    det_beta = _determinant(overlap_beta)
    det_total = det_alpha * det_beta
    sign = -1 if det_total < 0.0 else 1

    return DeterminantResult(
        alpha=det_alpha,
        beta=det_beta,
        total=det_total,
        sign=sign,
    )


def _write_nan_row(handle: Any, frame_id: int, time_fs: float) -> None:
    handle.write(
        f"{frame_id:5d}  {time_fs:10.3f}  "
        + "  ".join(["nan"] * 11)
        + "\n"
    )


def load_cdftbci_simple_config(config_path: Path) -> SimpleCDFTBCIConfig:
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent

    with open(config_path, "r") as handle:
        data = yaml.safe_load(handle)

    input_cfg = data.get("input", {})
    required_inputs = [
        "trajectory",
        "topology",
        "qm_atoms",
        "pccharges_template",
        "hsd_template",
    ]
    for key in required_inputs:
        if key not in input_cfg:
            raise ValueError(f"Missing required input: {key}")

    output_cfg = data.get("output", {})
    output_dir = base_dir / output_cfg.get("directory", "output")

    frames_cfg = data.get("frames", {})
    dftb_cfg = data.get("dftb", {})
    ci_cfg = data.get("cdftb_ci", {})
    scc_cfg = data.get("scc", {})
    phase_cfg = data.get("phase_tracking", {})
    odin_cfg = phase_cfg.get("odin", {})

    fragments: List[FragmentConfig] = []
    for frag_data in data.get("fragments", []):
        atom_range = frag_data["atom_range"]
        charge_sum_range = frag_data.get("charge_sum_range", parse_atom_range(atom_range))
        initial_charges = None
        if "initial_charges" in frag_data:
            initial_charges = base_dir / frag_data["initial_charges"]
        fragments.append(
            FragmentConfig(
                name=frag_data["name"],
                atom_range=atom_range,
                charge_sum_range=charge_sum_range,
                initial_charges=initial_charges,
            )
        )

    if len(fragments) != 2:
        raise ValueError("cdftbci_simple requires exactly two fragments")

    retry_cfg = scc_cfg.get("retry", None)
    retry_attempts: Optional[List[Dict[str, Any]]]
    if retry_cfg is None:
        retry_attempts = None
    else:
        if not isinstance(retry_cfg, list) or not retry_cfg:
            raise ValueError("scc.retry must be a non-empty list of mappings")
        retry_attempts = []
        for index, item in enumerate(retry_cfg):
            if not isinstance(item, dict):
                raise ValueError(f"scc.retry[{index}] must be a mapping")
            retry_attempts.append(
                {
                    "label": str(item.get("label", f"attempt{index}")),
                    "use_initial_charges": item.get("use_initial_charges", None),
                    "disable_constraint": bool(item.get("disable_constraint", False)),
                    "mixing_parameter": item.get("mixing_parameter", None),
                    "max_scc_iterations": item.get("max_scc_iterations", None),
                    "mixer": item.get("mixer", None),
                    "optimiser": item.get("optimiser", None),
                    "max_constr_iterations": item.get("max_constr_iterations", None),
                }
            )

    odin_executable = odin_cfg.get("executable")
    if not odin_executable:
        raise ValueError("phase_tracking.odin.executable is required")

    return SimpleCDFTBCIConfig(
        traj_path=base_dir / input_cfg["trajectory"],
        topology_path=base_dir / input_cfg["topology"],
        qm_atoms_file=base_dir / input_cfg["qm_atoms"],
        pccharges_template=base_dir / input_cfg["pccharges_template"],
        hsd_template=base_dir / input_cfg["hsd_template"],
        output_dir=output_dir,
        energy_file=output_dir / output_cfg.get("energy_file", "energies.dat"),
        charge_file=output_dir / output_cfg.get("charge_file", "charges.dat"),
        start_frame=int(frames_cfg.get("start", 0)),
        n_frames=frames_cfg.get("n_frames", None),
        t0_fs=float(frames_cfg.get("t0_fs", 0.0)),
        dt_fs=float(frames_cfg.get("dt_fs", 4.0)),
        dftb_library_path=dftb_cfg.get(
            "library_path",
            "/home/takahashi/opt/dftb+/lib/libdftbplus.so",
        ),
        num_threads=dftb_cfg.get("num_threads", None),
        timeout=int(dftb_cfg.get("timeout", 300)),
        fragments=fragments,
        ci_N_A=float(ci_cfg.get("N_A", 101.0)),
        ci_N_B=float(ci_cfg.get("N_B", 101.0)),
        ci_M_A=float(ci_cfg.get("M_A", 0.0)),
        ci_M_B=float(ci_cfg.get("M_B", 0.0)),
        work_directory=output_cfg.get("work_directory", "work"),
        output_file=output_dir / output_cfg.get("simple_file", "cdftbci_simple.dat"),
        use_previous_charges=bool(scc_cfg.get("use_previous_charges", False)),
        retry_attempts=retry_attempts,
        odin_executable=str(odin_executable),
        odin_sk_prefix=str(odin_cfg.get("sk_prefix", "")),
        odin_sk_separator=str(odin_cfg.get("sk_separator", "-")),
        odin_sk_suffix=str(odin_cfg.get("sk_suffix", ".skf")),
        odin_lmax={str(key): int(value) for key, value in odin_cfg.get("lmax", {}).items()},
        odin_work_subdir=str(odin_cfg.get("work_subdir", "odin_work")),
        odin_keep_files=bool(odin_cfg.get("keep_files", False)),
    )


def run_cdftbci_simple_analysis(config_path: Path) -> None:
    config = load_cdftbci_simple_config(config_path)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    qm_indices = load_qm_indices(config.qm_atoms_file)
    traj_info = md.load(str(config.traj_path), top=str(config.topology_path), frame=0)
    mm_indices = load_mm_indices(qm_indices, traj_info.n_atoms)
    mm_charges = load_pccharges(config.pccharges_template)
    atom_types = extract_atom_types_from_hsd(config.hsd_template)

    work_dir = config.output_dir / config.work_directory
    active_retry = config.retry_attempts if config.retry_attempts is not None else DEFAULT_RETRY_ATTEMPTS

    previous_coords_ang: Optional[np.ndarray] = None
    previous_state_a: Optional[GaugeState] = None
    previous_state_b: Optional[GaugeState] = None
    processed_count = 0

    with open(config.output_file, "w") as output:
        output.write(
            f"# {'frame':>5s}  {'time_fs':>10s}  {'eta_A':>6s}  {'eta_B':>6s}  {'gauge':>6s}  "
            f"{'D_A_alpha':>16s}  {'D_A_beta':>16s}  {'D_A':>16s}  "
            f"{'D_B_alpha':>16s}  {'D_B_beta':>16s}  {'D_B':>16s}  "
            f"{'S_AB':>16s}  {'J_lowdin_meV':>16s}\n"
        )

        for frame_id, _, qm_coords_bohr, mm_coords_ang in iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices, mm_indices
        ):
            if frame_id < config.start_frame:
                continue
            if config.n_frames is not None and processed_count >= config.n_frames:
                break

            qm_coords_ang = qm_coords_bohr / BOHR_PER_ANG
            time_fs = config.t0_fs + frame_id * config.dt_fs

            print(f"Frame {frame_id:05d} (t = {time_fs:.3f} fs)")
            setup_work_directory(
                work_dir,
                qm_coords_ang,
                mm_coords_ang,
                mm_charges,
                frame_id,
                time_fs,
                atom_types,
            )

            energies: List[float] = []
            success_all = True

            for frag in config.fragments:
                charges_file = work_dir / frag.name / "charges.dat"
                use_initial_charges = False

                if config.use_previous_charges:
                    if processed_count == 0:
                        if frag.initial_charges is not None and frag.initial_charges.exists():
                            use_initial_charges = True
                            charges_file.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy(frag.initial_charges, charges_file)
                    elif charges_file.exists():
                        use_initial_charges = True

                energy, _, success, retry_info, _ = run_cdftb_with_retry(
                    work_dir,
                    frag,
                    qm_coords_bohr,
                    config.hsd_template,
                    config.dftb_library_path,
                    config.num_threads,
                    config.timeout,
                    use_initial_charges,
                    True,
                    retry_attempts=active_retry,
                )

                if not success:
                    success_all = False
                    print(f"  {frag.name}: ERROR - {retry_info}")
                    break

                print(f"  {frag.name}: E = {energy:.10f} a.u.")
                energies.append(float(energy))

                run_dftb_in_subprocess(
                    qm_coords_bohr,
                    work_dir / frag.name,
                    write_hs=True,
                    dftb_library_path=config.dftb_library_path,
                    num_threads=config.num_threads,
                    timeout=config.timeout,
                )

            if not success_all:
                _write_nan_row(output, frame_id, time_fs)
                output.flush()
                processed_count += 1
                continue

            ham, _, _, error, orb_a, orb_b, _ = compute_cdftbci_for_frame(
                work_dir,
                config.fragments[0].name,
                config.fragments[1].name,
                energies[0],
                energies[1],
                config.ci_N_A,
                config.ci_N_B,
                config.ci_M_A,
                config.ci_M_B,
            )

            if error or ham is None or orb_a is None or orb_b is None:
                print(f"  CDFTB-CI: ERROR - {error}")
                _write_nan_row(output, frame_id, time_fs)
                output.flush()
                processed_count += 1
                continue

            det_a = DeterminantResult(np.nan, np.nan, np.nan, 1)
            det_b = DeterminantResult(np.nan, np.nan, np.nan, 1)
            eta_a = 1
            eta_b = 1

            if (
                previous_coords_ang is not None
                and previous_state_a is not None
                and previous_state_b is not None
            ):
                try:
                    odin_work_dir = None
                    if config.odin_work_subdir:
                        odin_work_dir = config.output_dir / config.odin_work_subdir
                        if config.odin_keep_files:
                            odin_work_dir = odin_work_dir / f"frame_{frame_id:05d}"

                    s_ao_cross = compute_cross_overlap_odin(
                        previous_coords_ang,
                        qm_coords_ang,
                        atom_types,
                        config.odin_lmax,
                        config.odin_sk_prefix,
                        config.odin_sk_separator,
                        config.odin_sk_suffix,
                        config.odin_executable,
                        remove_translation_rotation=False,
                        work_dir=odin_work_dir,
                        keep_files=config.odin_keep_files,
                    )
                    det_a = _compute_state_determinant(previous_state_a, orb_a, s_ao_cross)
                    det_b = _compute_state_determinant(previous_state_b, orb_b, s_ao_cross)
                    eta_a = det_a.sign
                    eta_b = det_b.sign
                except Exception as exc:
                    print(f"  ODIN overlap: ERROR - {exc}")

            gauge = eta_a * eta_b
            if previous_state_a is None or previous_state_b is None:
                j_lowdin_raw_meV = 27211.386 * compute_transfer_integral_unrestricted(
                    ham.H,
                    ham.S,
                    method="lowdin",
                )
                if np.isfinite(j_lowdin_raw_meV) and j_lowdin_raw_meV < 0.0:
                    eta_a = -1
                    gauge = -1

            ham.H_AB = gauge * ham.H_AB
            ham.S_AB = gauge * ham.S_AB
            ham.H = np.array([[ham.E_A, ham.H_AB], [ham.H_AB, ham.E_B]])
            ham.S = np.array([[1.0, ham.S_AB], [ham.S_AB, 1.0]])

            j_lowdin_meV = 27211.386 * compute_transfer_integral_unrestricted(
                ham.H,
                ham.S,
                method="lowdin",
            )

            output.write(
                f"{frame_id:5d}  {time_fs:10.3f}  "
                f"{eta_a:6d}  {eta_b:6d}  {gauge:6d}  "
                f"{det_a.alpha:16.10f}  {det_a.beta:16.10f}  {det_a.total:16.10f}  "
                f"{det_b.alpha:16.10f}  {det_b.beta:16.10f}  {det_b.total:16.10f}  "
                f"{ham.S_AB:16.10f}  {j_lowdin_meV:16.10f}\n"
            )
            output.flush()

            previous_coords_ang = qm_coords_ang.copy()
            previous_state_a = _apply_state_sign(orb_a, eta_a)
            previous_state_b = _apply_state_sign(orb_b, eta_b)
            processed_count += 1

            print(
                f"  CDFTB-CI: D_A={det_a.total:+.6f}, D_B={det_b.total:+.6f}, "
                f"S_AB={ham.S_AB:+.6f}, J_lowdin={j_lowdin_meV:+.4f} meV"
            )

    print(f"Saved simple CDFTB-CI results to {config.output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Minimal CDFTB-CI analysis with ODIN-based phase correction")
    parser.add_argument("config", type=Path, help="Path to YAML configuration file")
    args = parser.parse_args()
    run_cdftbci_simple_analysis(args.config)


if __name__ == "__main__":
    main()