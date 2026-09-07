"""Run NWChem constrained-DFT (CDFT) calculations along a trajectory.

This module generates NWChem inputs with charge and/or spin constraints,
executes NWChem in a subprocess for each trajectory frame, and extracts the
total DFT energy and constrained populations.

The CDFT directives follow NWChem's ``dft`` block syntax::

    dft
      xc b3lyp
      mult 1
      odft
      convergence nolevelshifting
      cdft 1 36 charge 1.0 pop lowdin
      cdft 1 36 spin 1.0 pop lowdin
    end

Point charges for QM/MM calculations can be appended to the geometry block as
``bq`` records.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import yaml

ANG_PER_NM = 10.0


def parse_atom_range(atom_range: str) -> List[int]:
    """Parse a NWChem-style, 1-indexed atom range into Python ``[start, end)``.

    Examples
    --------
    ``"1:36"`` -> ``[0, 36]``
    ``"37:72"`` -> ``[36, 72]``
    """
    parts = str(atom_range).split(":")
    if len(parts) != 2:
        raise ValueError(f"Invalid atom range {atom_range!r}; expected 'start:end'")
    try:
        first = int(parts[0])
        last = int(parts[1])
    except ValueError as exc:
        raise ValueError(f"Invalid atom range {atom_range!r}") from exc
    if first < 1 or last < first:
        raise ValueError(f"Invalid atom range {atom_range!r}; expected 1 <= start <= end")
    return [first - 1, last]


@dataclass
class NWChemConstraint:
    """A single NWChem CDFT constraint."""

    kind: str = "charge"  # "charge" or "spin"
    atoms: str = ""       # 1-indexed range, e.g. "1:36"
    value: float = 0.0
    population: str = "lowdin"

    def __post_init__(self) -> None:
        kind = self.kind.lower()
        if kind not in {"charge", "spin"}:
            raise ValueError(f"Unsupported CDFT constraint kind: {self.kind!r}")
        self.kind = kind
        population = self.population.lower()
        if population not in {"lowdin", "becke", "mulliken"}:
            raise ValueError(f"Unsupported CDFT population scheme: {self.population!r}")
        self.population = population


@dataclass
class NWChemCDFTConfig:
    """Configuration for NWChem constrained-DFT trajectory analysis."""

    traj_path: Path
    topology_path: Path
    qm_atoms_file: Path
    basis: str = "6-31G*"
    xc: str = "b3lyp"
    charge: int = 0
    multiplicity: int = 1

    pointcharges_template: Optional[Path] = None

    output_dir: Path = Path("output_nwchem")
    energy_file: Path = Path("output_nwchem/energies.dat")
    population_file: Path = Path("output_nwchem/populations.dat")
    work_directory: str = "work"

    start_frame: int = 0
    n_frames: Optional[int] = None
    t0_fs: float = 0.0
    dt_fs: float = 1.0

    nwchem_binary: str = "nwchem"
    scratch_dir: Optional[Path] = None
    timeout: int = 3600

    constraints: List[NWChemConstraint] = field(default_factory=list)


def build_cdft_directive(constraint: NWChemConstraint) -> str:
    """Render one NWChem ``cdft`` directive."""
    start, end = parse_atom_range(constraint.atoms)
    # parse_atom_range returns Python [start, end); NWChem uses inclusive 1-based
    # atom numbers, so start+1 .. end.
    value = float(constraint.value)
    if value.is_integer():
        value_text = f"{value:.1f}"
    else:
        value_text = f"{value:g}"
    return (
        f"cdft {start + 1} {end} {constraint.kind} "
        f"{value_text} pop {constraint.population}"
    )


def render_nwchem_input(
    atom_types: Sequence[str],
    coords_ang: np.ndarray | Sequence[Sequence[float]],
    constraints: Sequence[NWChemConstraint],
    basis: str = "6-31G*",
    xc: str = "b3lyp",
    charge: int = 0,
    multiplicity: int = 1,
    point_charges: Optional[np.ndarray | Sequence[Sequence[float]]] = None,
    print_mo_vectors: bool = False,
) -> str:
    """Render a complete NWChem CDFT input deck.

    Parameters
    ----------
    atom_types
        Element symbols, one per QM atom.
    coords_ang
        QM coordinates in Angstrom with shape ``(n_atoms, 3)``.
    constraints
        CDFT constraints to place inside the ``dft`` block.
    point_charges
        Optional ``(n, 4)`` array of ``x y z charge`` records in Angstrom.
    """
    coords = np.asarray(coords_ang, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords_ang must have shape (n_atoms, 3)")
    if len(atom_types) != coords.shape[0]:
        raise ValueError("atom_types and coords_ang have different numbers of atoms")

    lines: List[str] = []
    lines.append("start molecule")
    lines.append('title "q4mdkit NWChem CDFT"')
    lines.append(f"charge {int(charge)}")
    lines.append("geometry noautoz nocenter noautosym")
    for atom, xyz in zip(atom_types, coords):
        lines.append(
            f"{atom:2s} {xyz[0]:15.10f} {xyz[1]:15.10f} {xyz[2]:15.10f}"
        )
    if point_charges is not None:
        pc = np.asarray(point_charges, dtype=float)
        if pc.ndim != 2 or pc.shape[1] != 4:
            raise ValueError("point_charges must have shape (n, 4)")
        for x, y, z, q in pc:
            lines.append(f"{'bq':3s} {x:15.10f} {y:15.10f} {z:15.10f} {q:15.10f}")
    lines.append("end")
    lines.append("")
    lines.append("basis")
    lines.append(f"  * library {basis}")
    lines.append("end")
    lines.append("")
    if print_mo_vectors:
        lines.append("set movecs:tanalyze 0.0")
        lines.append("")
    lines.append("dft")
    lines.append(f"  xc {xc}")
    lines.append(f"  mult {int(multiplicity)}")
    lines.append("  odft")
    lines.append("  convergence nolevelshifting")
    for constraint in constraints:
        lines.append(f"  {build_cdft_directive(constraint)}")
    if print_mo_vectors:
        lines.append('  print "final vectors analysis" "final vectors"')
    lines.append("end")
    lines.append("")
    lines.append("task dft energy")
    return "\n".join(lines) + "\n"


def parse_nwchem_energy(output_text: str) -> float:
    """Extract the final total DFT energy (Hartree) from NWChem output."""
    match = re.search(
        r"Total\s+DFT\s+energy\s*=\s*([-+]?\d+\.\d+)",
        output_text,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError("Could not find 'Total DFT energy' in NWChem output")
    return float(match.group(1))


def parse_nwchem_populations(output_text: str) -> Optional[List[float]]:
    """Extract the CDFT constraint multipliers printed by NWChem.

    NWChem reports the final Lagrange multiplier for each constraint under a
    ``CDFT final multipliers`` block. The values are returned in constraint
    order. Older NWChem versions that print ``constraint value(s)`` are still
    accepted as a fallback. Returns ``None`` when no block can be found.
    """
    final = re.search(
        r"CDFT\s+final\s+multipliers\s*\n(.*?)(?=\n\s*\n|\Z)",
        output_text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if final:
        values = [float(v) for v in re.findall(r"[-+]?\d+\.\d+", final.group(1))]
        if values:
            return values

    patterns = [
        r"constraint\s+(?:value|values?)\s*[=:]\s*([-+]?\d+\.\d+)",
        r"CDFT\s+constrained\s+(?:value|values?)\s*[=:]\s*([-+]?\d+\.\d+)",
    ]
    for pattern in patterns:
        values = [float(v) for v in re.findall(pattern, output_text, flags=re.IGNORECASE)]
        if values:
            return values
    return None


def _resolve(path_text: str, base_dir: Path) -> Path:
    return (base_dir / path_text).resolve()


def load_config(config_path: str | Path) -> NWChemCDFTConfig:
    """Load a YAML configuration for NWChem CDFT trajectory analysis."""
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent

    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    input_cfg = data.get("input", {})
    for key in ("trajectory", "topology", "qm_atoms"):
        if key not in input_cfg:
            raise ValueError(f"Missing required input field: {key}")

    output_cfg = data.get("output", {})
    output_dir = _resolve(output_cfg.get("directory", "output_nwchem"), base_dir)
    energy_file = _resolve(output_cfg.get("energy_file", "energies.dat"), output_dir)
    population_file = _resolve(
        output_cfg.get("population_file", "populations.dat"), output_dir
    )
    work_directory = output_cfg.get("work_directory", "work")

    frames_cfg = data.get("frames", {})

    nwchem_cfg = data.get("nwchem", {})

    pointcharges_template = None
    if "pointcharges_template" in input_cfg:
        pointcharges_template = _resolve(input_cfg["pointcharges_template"], base_dir)

    constraints: List[NWChemConstraint] = []
    for item in data.get("constraints", []):
        constraints.append(
            NWChemConstraint(
                kind=str(item.get("kind", "charge")),
                atoms=str(item["atoms"]),
                value=float(item.get("value", 0.0)),
                population=str(item.get("population", "lowdin")),
            )
        )

    scratch_dir = None
    if "scratch_dir" in nwchem_cfg:
        scratch_dir = _resolve(nwchem_cfg["scratch_dir"], base_dir)

    return NWChemCDFTConfig(
        traj_path=_resolve(input_cfg["trajectory"], base_dir),
        topology_path=_resolve(input_cfg["topology"], base_dir),
        qm_atoms_file=_resolve(input_cfg["qm_atoms"], base_dir),
        basis=str(input_cfg.get("basis", "6-31G*")),
        xc=str(input_cfg.get("xc", "b3lyp")),
        charge=int(input_cfg.get("charge", 0)),
        multiplicity=int(input_cfg.get("multiplicity", 1)),
        pointcharges_template=pointcharges_template,
        output_dir=output_dir,
        energy_file=energy_file,
        population_file=population_file,
        work_directory=work_directory,
        start_frame=int(frames_cfg.get("start", 0)),
        n_frames=frames_cfg.get("n_frames"),
        t0_fs=float(frames_cfg.get("t0_fs", 0.0)),
        dt_fs=float(frames_cfg.get("dt_fs", 1.0)),
        nwchem_binary=str(nwchem_cfg.get("binary", "nwchem")),
        scratch_dir=scratch_dir,
        timeout=int(nwchem_cfg.get("timeout", 3600)),
        constraints=constraints,
    )


def _load_qm_indices(path: Path) -> np.ndarray:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"QM atom file not found: {path}")
    tokens = path.read_text().split()
    if not tokens:
        raise ValueError(f"QM atom file is empty: {path}")
    return np.array(sorted({int(token) for token in tokens}), dtype=int)


def _iter_qm_coordinates(
    traj_path: Path,
    topology_path: Path,
    qm_indices: np.ndarray,
    chunk_size: int = 10,
):
    """Yield ``(frame_id, time_fs, qm_coords_ang)`` using mdtraj lazily."""
    try:
        import mdtraj as md
    except ImportError as exc:
        raise ImportError("mdtraj is required for trajectory iteration") from exc

    frame_counter = 0
    for chunk in md.iterload(str(traj_path), top=str(topology_path), chunk=chunk_size):
        times = chunk.time if chunk.time is not None else [None] * chunk.n_frames
        qm_coords = chunk.xyz[:, qm_indices, :] * ANG_PER_NM
        for local_idx in range(chunk.n_frames):
            yield (
                frame_counter + local_idx,
                times[local_idx],
                np.asarray(qm_coords[local_idx], dtype=float),
            )
        frame_counter += chunk.n_frames


def _load_point_charges(path: Path) -> Optional[np.ndarray]:
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Point-charge file not found: {path}")
    return np.loadtxt(path)


def _extract_atom_types_from_topology(
    topology_path: Path,
    qm_indices: np.ndarray,
) -> List[str]:
    """Extract element symbols for QM atoms from a topology."""
    try:
        import mdtraj as md
    except ImportError as exc:
        raise ImportError("mdtraj is required to read topology atom types") from exc

    top = md.load_topology(str(topology_path))
    table, bonds = top.to_dataframe()
    if "element" not in table.columns:
        raise ValueError("Topology does not contain element information")
    symbols = table["element"].tolist()
    return [str(symbols[int(i)]) for i in qm_indices]


def _run_nwchem(
    work_dir: Path,
    config: NWChemCDFTConfig,
) -> subprocess.CompletedProcess[str]:
    cmd = [config.nwchem_binary, "nwchem.inp"]
    if config.scratch_dir is not None:
        cmd.append(str(config.scratch_dir))
    return subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=config.timeout,
        check=False,
    )


def run_nwchem_cdft_analysis(config_path: str | Path) -> None:
    """Run NWChem CDFT calculations for selected trajectory frames."""
    config = load_config(config_path)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    qm_indices = _load_qm_indices(config.qm_atoms_file)
    atom_types = _extract_atom_types_from_topology(config.topology_path, qm_indices)
    point_charges = _load_point_charges(config.pointcharges_template)

    work_dir = config.output_dir / config.work_directory
    energy_path = config.energy_file
    population_path = config.population_file

    with energy_path.open("w", encoding="utf-8") as energy_handle, population_path.open(
        "w", encoding="utf-8"
    ) as population_handle:
        energy_handle.write("# Frame  Time(fs)  TotalDFTEnergy(a.u.)\n")
        population_handle.write(
            "# Frame  Time(fs)  "
            + "  ".join(f"{c.kind}_{i + 1}" for i, c in enumerate(config.constraints))
            + "\n"
        )

        processed = 0
        for frame_id, traj_time_fs, qm_coords_ang in _iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices
        ):
            if frame_id < config.start_frame:
                continue
            if config.n_frames is not None and processed >= config.n_frames:
                break

            time_fs = config.t0_fs + frame_id * config.dt_fs
            work_dir.mkdir(parents=True, exist_ok=True)
            input_text = render_nwchem_input(
                atom_types=atom_types,
                coords_ang=qm_coords_ang,
                constraints=config.constraints,
                basis=config.basis,
                xc=config.xc,
                charge=config.charge,
                multiplicity=config.multiplicity,
                point_charges=point_charges,
            )
            (work_dir / "nwchem.inp").write_text(input_text, encoding="utf-8")

            completed = _run_nwchem(work_dir, config)
            output_text = completed.stdout + "\n" + completed.stderr

            if completed.returncode == 0:
                try:
                    energy = parse_nwchem_energy(output_text)
                except ValueError:
                    energy = float("nan")
                populations = parse_nwchem_populations(output_text)
            else:
                energy = float("nan")
                populations = None

            energy_handle.write(f"{frame_id:5d}  {time_fs:12.3f}  {energy:18.10f}\n")
            energy_handle.flush()

            if populations is None:
                populations = [float("nan")] * len(config.constraints)
            if len(populations) != len(config.constraints):
                populations = populations[: len(config.constraints)]
                populations.extend(
                    [float("nan")] * (len(config.constraints) - len(populations))
                )
            population_values = "  ".join(f"{value:12.8f}" for value in populations)
            population_handle.write(
                f"{frame_id:5d}  {time_fs:12.3f}  {population_values}\n"
            )
            population_handle.flush()
            processed += 1


def main() -> None:
    """Command-line entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Run NWChem constrained-DFT analysis")
    parser.add_argument("config", help="Path to YAML configuration file")
    args = parser.parse_args()
    run_nwchem_cdft_analysis(args.config)


if __name__ == "__main__":
    main()
