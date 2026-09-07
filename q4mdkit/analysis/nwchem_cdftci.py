"""NWChem CDFT-CI backend.

Extracts the pieces needed for the shared CDFT-CI routines
(``cdft_ci_core``) from NWChem text output and drives the two diabatic
constrained-DFT states per trajectory frame.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple

import numpy as np
import yaml

from .cdft_ci_core import UnrestrictedOrbitalData
from .cdft_ci_core import (
    UnrestrictedCDFTCIHamiltonian,
    build_cdft_ci_hamiltonian_unrestricted,
    build_lowdin_weight_matrix,
    build_mulliken_weight_matrix,
    compute_coupling_element_unrestricted,
    compute_transfer_integral_unrestricted,
)
from .nwchem_cdft import (
    NWChemCDFTConfig,
    NWChemConstraint,
    _extract_atom_types_from_topology,
    _iter_qm_coordinates,
    _load_qm_indices,
    _run_nwchem,
    load_config,
    parse_nwchem_energy,
    parse_nwchem_populations,
    render_nwchem_input,
)


_VECTOR_RE = re.compile(
    r"^\s*Vector\s+(\d+)\s+Occ\s*=\s*([-+]?\d+(?:\.\d+)?(?:[DdEe][+-]?\d+)?)",
    re.IGNORECASE,
)
_BFN_COEFF_RE = re.compile(r"(\d+)\s+([-+]?\d+\.\d+)")
_ATOM_BFN_RE = re.compile(r"(\d+)\s+([-+]?\d+\.\d+)\s+(\d+)")
_EVECS_HEADER_RE = re.compile(
    r"global array:\s*(alpha|beta)\s+evecs\[1:(\d+),1:(\d+)\]",
    re.IGNORECASE,
)

_ATOMIC_NUMBERS = {
    "H": 1, "He": 2,
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Br": 35, "I": 53,
}


def _to_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _extract_spin_section(output_text: str, spin: str) -> List[str]:
    header = f"DFT Final {spin} Molecular Orbital Analysis"
    lines = output_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if header in line:
            start = i
            break
    if start is None:
        raise ValueError(f"Could not find {header!r} in NWChem output")

    end = len(lines)
    for i in range(start + 1, len(lines)):
        if "DFT Final" in lines[i] or "alpha - beta orbital overlaps" in lines[i]:
            end = i
            break
    return lines[start:end]


def _parse_spin_section(
    output_text: str, spin: str
) -> Tuple[np.ndarray, np.ndarray]:
    section = _extract_spin_section(output_text, spin)

    n_mo = 0
    occ: List[float] = []
    entries: List[Tuple[int, int, float]] = []  # (mo_idx, bfn_idx, coeff)

    current_mo = -1
    for line in section:
        vec_match = _VECTOR_RE.search(line)
        if vec_match:
            current_mo = int(vec_match.group(1)) - 1
            n_mo = max(n_mo, int(vec_match.group(1)))
            occ.append(_to_float(vec_match.group(2)))
            continue

        if current_mo < 0:
            continue
        if "Bfn." in line or "MO Center" in line or line.strip().startswith("---"):
            continue
        if not line.strip():
            continue

        for bfn_text, coeff_text in _BFN_COEFF_RE.findall(line):
            entries.append((current_mo, int(bfn_text) - 1, float(coeff_text)))

    if n_mo == 0:
        raise ValueError(f"No molecular orbitals found for {spin} spin")

    n_bf = max(e[1] for e in entries) + 1
    C = np.zeros((n_bf, n_mo), dtype=np.float64)
    for mo_idx, bfn_idx, coeff in entries:
        C[bfn_idx, mo_idx] = coeff

    occ_array = np.array(occ, dtype=np.float64)
    if occ_array.shape[0] < n_mo:
        occ_array = np.pad(occ_array, (0, n_mo - occ_array.shape[0]))
    return C, occ_array


def _parse_evecs_matrix(output_text: str, spin: str) -> np.ndarray:
    """Parse the true normalized MO coefficient matrix from ``Final MO vectors``."""
    header = re.search(
        rf"global array:\s*{spin}\s+evecs\[1:(\d+),1:(\d+)\]",
        output_text,
        flags=re.IGNORECASE,
    )
    if not header:
        raise ValueError(f"Final MO vectors for {spin} not found in NWChem output")

    n_bf = int(header.group(1))
    n_mo = int(header.group(2))
    lines = output_text.splitlines()

    start = None
    for i, line in enumerate(lines):
        if re.search(rf"global array:\s*{spin}\s+evecs", line, flags=re.IGNORECASE):
            start = i
            break
    if start is None:
        raise ValueError(f"Final MO vectors for {spin} not found")

    C = np.zeros((n_bf, n_mo), dtype=np.float64)
    i = start + 1
    for _block in range((n_mo + 5) // 6):
        # Locate the column-index header line.
        while i < len(lines):
            tokens = lines[i].split()
            if tokens and all(t.replace("-", "").isdigit() for t in tokens):
                col_indices = [int(t) - 1 for t in tokens]
                i += 1
                break
            i += 1
        # Skip the separator line and any blank lines.
        while i < len(lines) and (
            not lines[i].strip() or lines[i].strip().startswith("-")
        ):
            i += 1
        for _row in range(n_bf):
            while i < len(lines) and not lines[i].strip():
                i += 1
            tokens = lines[i].split()
            row = int(tokens[0]) - 1
            values = [float(t) for t in tokens[1:]]
            for col, value in zip(col_indices, values):
                C[row, col] = value
            i += 1
    return C


def parse_nwchem_mo_coefficients(
    output_text: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(C_alpha, C_beta, occ_alpha, occ_beta)`` from NWChem output.

    The true normalized MO coefficients are read from the ``Final MO vectors``
    blocks. Occupations are read from the ``Molecular Orbital Analysis``
    sections.
    """
    C_alpha = _parse_evecs_matrix(output_text, "alpha")
    C_beta = _parse_evecs_matrix(output_text, "beta")
    _, occ_alpha = _parse_spin_section(output_text, "Alpha")
    _, occ_beta = _parse_spin_section(output_text, "Beta")
    return C_alpha, C_beta, occ_alpha, occ_beta


def reconstruct_ao_overlap(C: np.ndarray) -> np.ndarray:
    """Reconstruct the AO overlap from orthonormal MO coefficients.

    With ``C^T S C = I`` and a square, full-rank ``C``,
    ``S = (C C^T)^{-1}``.
    """
    gram = C @ C.T
    try:
        return np.linalg.inv(gram)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(gram)


def build_nwchem_orbital_data(
    output_text: str,
) -> UnrestrictedOrbitalData:
    """Build ``UnrestrictedOrbitalData`` from an NWChem CDFT output text."""
    C_alpha, C_beta, occ_alpha, occ_beta = parse_nwchem_mo_coefficients(output_text)
    n_alpha = int(np.sum(occ_alpha > 0.5))
    n_beta = int(np.sum(occ_beta > 0.5))
    return UnrestrictedOrbitalData(
        C_alpha=C_alpha,
        C_beta=C_beta,
        occ_alpha=occ_alpha,
        occ_beta=occ_beta,
        n_alpha=n_alpha,
        n_beta=n_beta,
    )


def parse_nwchem_atom_orbital_map(output_text: str) -> dict[int, list[int]]:
    """Return a 1-based atom index -> 0-based AO indices mapping."""
    section = _extract_spin_section(output_text, "Alpha")
    atom_to_orbitals: dict[int, list[int]] = {}
    for line in section:
        if "Bfn." in line or "MO Center" in line or line.strip().startswith("---"):
            continue
        for bfn_text, _coeff, atom_text in _ATOM_BFN_RE.findall(line):
            bfn = int(bfn_text) - 1
            atom = int(atom_text)
            atom_to_orbitals.setdefault(atom, [])
            if bfn not in atom_to_orbitals[atom]:
                atom_to_orbitals[atom].append(bfn)
    return atom_to_orbitals


def compute_nwchem_cdftci_from_outputs(
    output_A: str,
    output_B: str,
    fragment_A_atoms: list[int],
    fragment_B_atoms: list[int],
    N_A: float,
    N_B: float,
    M_A: float = 0.0,
    M_B: float = 0.0,
    Z_A: float = 0.0,
    Z_B: float = 0.0,
    weight_scheme: str = "lowdin",
) -> UnrestrictedCDFTCIHamiltonian:
    """Compute the CDFT-CI Hamiltonian from two NWChem CDFT output texts."""
    orb_A = build_nwchem_orbital_data(output_A)
    orb_B = build_nwchem_orbital_data(output_B)

    S_AO = reconstruct_ao_overlap(orb_A.C_alpha)
    atom_to_orbitals = parse_nwchem_atom_orbital_map(output_A)
    n_atoms = max(atom_to_orbitals)

    if weight_scheme == "lowdin":
        build_weight = build_lowdin_weight_matrix
    elif weight_scheme == "mulliken":
        build_weight = build_mulliken_weight_matrix
    else:
        raise ValueError(f"Unknown weight scheme: {weight_scheme!r}")

    w_A = build_weight(S_AO, fragment_A_atoms, n_atoms, atom_to_orbitals)
    w_B = build_weight(S_AO, fragment_B_atoms, n_atoms, atom_to_orbitals)

    E_A = parse_nwchem_energy(output_A)
    E_B = parse_nwchem_energy(output_B)

    mult_A = parse_nwchem_populations(output_A) or [0.0]
    mult_B = parse_nwchem_populations(output_B) or [0.0]
    V_A = mult_A[0]
    V_B = mult_B[0]
    V_M_A = mult_A[1] if len(mult_A) > 1 else 0.0
    V_M_B = mult_B[1] if len(mult_B) > 1 else 0.0

    ham = build_cdft_ci_hamiltonian_unrestricted(
        orb_A,
        orb_B,
        S_AO,
        w_A,
        w_B,
        E_A,
        E_B,
        V_A,
        V_B,
        N_A,
        N_B,
        V_M_A,
        V_M_B,
        M_A,
        M_B,
    )

    # NWChem's `cdft ... charge` constrains the net charge, while the Mulliken
    # weight matrix above measures electron population. Convert the charge
    # weight overlaps: <1|Z - n|2> = Z * S_AB - <1|n|2>.
    if Z_A != 0.0 or Z_B != 0.0:
        S_AB = ham.S_AB
        W_BA_charge = Z_A * S_AB - ham.W_BA
        W_AB_charge = Z_B * S_AB - ham.W_AB
        H_AB = compute_coupling_element_unrestricted(
            E_A, E_B, V_A, V_B, N_A, N_B, S_AB, W_BA_charge, W_AB_charge,
            V_M_A, V_M_B, M_A, M_B, ham.W_M_BA, ham.W_M_AB,
        )
        ham.W_BA = W_BA_charge
        ham.W_AB = W_AB_charge
        ham.H_AB = H_AB
        ham.H = np.array([[E_A, H_AB], [H_AB, E_B]])
        ham.S = np.array([[1.0, S_AB], [S_AB, 1.0]])

    return ham


@dataclass
class NWChemCDFTCIState:
    """A single diabatic CDFT-CI state."""

    name: str
    atoms: str  # NWChem 1-based inclusive range, e.g. "1:36"
    charge: float  # target charge population N
    spin: float = 0.0  # target spin population M (0.0 -> no spin constraint)
    population: str = "lowdin"

    def build_constraints(self) -> List[NWChemConstraint]:
        constraints = [
            NWChemConstraint(
                kind="charge",
                atoms=self.atoms,
                value=self.charge,
                population=self.population,
            )
        ]
        if self.spin != 0.0:
            constraints.append(
                NWChemConstraint(
                    kind="spin",
                    atoms=self.atoms,
                    value=self.spin,
                    population=self.population,
                )
            )
        return constraints


@dataclass
class NWChemCDFTCIConfig(NWChemCDFTConfig):
    """Configuration for a two-state NWChem CDFT-CI trajectory analysis."""

    state_A: NWChemCDFTCIState = field(
        default_factory=lambda: NWChemCDFTCIState("state_a", "1:36", 1.0)
    )
    state_B: NWChemCDFTCIState = field(
        default_factory=lambda: NWChemCDFTCIState("state_b", "37:72", -1.0)
    )
    coupling_file: Path = Path("output_nwchem/cdftci.dat")


def load_cdftci_config(config_path: str | Path) -> NWChemCDFTCIConfig:
    """Load a YAML configuration for NWChem CDFT-CI trajectory analysis."""
    base = load_config(config_path)
    config_path = Path(config_path).resolve()
    base_dir = config_path.parent

    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    ci_cfg = data.get("cdftci", {})

    def _state(key: str, default_atoms: str, default_charge: float) -> NWChemCDFTCIState:
        item = ci_cfg.get(key, {})
        return NWChemCDFTCIState(
            name=str(item.get("name", key)),
            atoms=str(item.get("atoms", default_atoms)),
            charge=float(item.get("charge", default_charge)),
            spin=float(item.get("spin", 0.0)),
            population=str(item.get("population", "mulliken")),
        )

    coupling_file = base.output_dir / ci_cfg.get("coupling_file", "cdftci.dat")

    return NWChemCDFTCIConfig(
        traj_path=base.traj_path,
        topology_path=base.topology_path,
        qm_atoms_file=base.qm_atoms_file,
        basis=base.basis,
        xc=base.xc,
        charge=base.charge,
        multiplicity=base.multiplicity,
        pointcharges_template=base.pointcharges_template,
        output_dir=base.output_dir,
        energy_file=base.energy_file,
        population_file=base.population_file,
        work_directory=base.work_directory,
        start_frame=base.start_frame,
        n_frames=base.n_frames,
        t0_fs=base.t0_fs,
        dt_fs=base.dt_fs,
        nwchem_binary=base.nwchem_binary,
        scratch_dir=base.scratch_dir,
        timeout=base.timeout,
        state_A=_state("state_a", "1:36", 1.0),
        state_B=_state("state_b", "37:72", -1.0),
        coupling_file=coupling_file,
    )


def _run_state(
    work_dir: Path,
    config: NWChemCDFTCIConfig,
    atom_types: List[str],
    coords_ang: np.ndarray,
    state: NWChemCDFTCIState,
) -> str:
    work_dir.mkdir(parents=True, exist_ok=True)
    text = render_nwchem_input(
        atom_types=atom_types,
        coords_ang=coords_ang,
        constraints=state.build_constraints(),
        basis=config.basis,
        xc=config.xc,
        charge=config.charge,
        multiplicity=config.multiplicity,
        point_charges=None,
        print_mo_vectors=True,
    )
    (work_dir / "nwchem.inp").write_text(text, encoding="utf-8")
    completed = _run_nwchem(work_dir, config)
    output_text = completed.stdout + "\n" + completed.stderr
    (work_dir / "nwchem.out").write_text(output_text, encoding="utf-8")
    return output_text


def run_nwchem_cdftci_analysis(config_path: str | Path) -> None:
    """Run two-state NWChem CDFT-CI calculations along a trajectory."""
    config = load_cdftci_config(config_path)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    qm_indices = _load_qm_indices(config.qm_atoms_file)
    atom_types = _extract_atom_types_from_topology(config.topology_path, qm_indices)
    state_A = config.state_A
    state_B = config.state_B

    fragment_A_atoms = _fragment_atom_list(state_A.atoms)
    fragment_B_atoms = _fragment_atom_list(state_B.atoms)
    Z_A = _fragment_nuclear_charge(atom_types, fragment_A_atoms)
    Z_B = _fragment_nuclear_charge(atom_types, fragment_B_atoms)

    coupling_path = config.coupling_file
    coupling_path.parent.mkdir(parents=True, exist_ok=True)
    with coupling_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Frame  Time(fs)  E_A(Ha)  E_B(Ha)  S_AB  H_AB(Ha)  J_direct(Ha)  "
            "J_lowdin(Ha)\n"
        )

        processed = 0
        frames = _iter_qm_coordinates(
            config.traj_path, config.topology_path, qm_indices
        )
        for frame_id, _traj_time, qm_coords_ang in frames:
            if frame_id < config.start_frame:
                continue
            if config.n_frames is not None and processed >= config.n_frames:
                break

            time_fs = config.t0_fs + frame_id * config.dt_fs
            work_dir = config.output_dir / config.work_directory
            print(
                f"[CDFT-CI] frame {frame_id}: running {state_A.name} then {state_B.name}",
                flush=True,
            )
            out_A = _run_state(
                work_dir / state_A.name,
                config,
                atom_types,
                qm_coords_ang,
                state_A,
            )
            print(f"[CDFT-CI] frame {frame_id}: {state_A.name} done", flush=True)
            out_B = _run_state(
                work_dir / state_B.name,
                config,
                atom_types,
                qm_coords_ang,
                state_B,
            )
            print(f"[CDFT-CI] frame {frame_id}: {state_B.name} done", flush=True)

            try:
                ham = compute_nwchem_cdftci_from_outputs(
                    out_A,
                    out_B,
                    fragment_A_atoms,
                    fragment_B_atoms,
                    N_A=state_A.charge,
                    N_B=state_B.charge,
                    M_A=state_A.spin,
                    M_B=state_B.spin,
                    Z_A=Z_A,
                    Z_B=Z_B,
                    weight_scheme=state_A.population,
                )
                J_direct = compute_transfer_integral_unrestricted(
                    ham.H, ham.S, method="direct"
                )
                J_lowdin = compute_transfer_integral_unrestricted(
                    ham.H, ham.S, method="lowdin"
                )
                handle.write(
                    f"{frame_id:6d}  {time_fs:8.2f}  {ham.E_A:14.10f}  "
                    f"{ham.E_B:14.10f}  {ham.S_AB:12.8f}  {ham.H_AB:14.10f}  "
                    f"{J_direct:14.10f}  {J_lowdin:14.10f}\n"
                )
                handle.flush()
            except Exception as exc:  # keep going; write NaNs
                handle.write(
                    f"{frame_id:6d}  {time_fs:8.2f}  nan  nan  nan  nan  nan  nan"
                    f"  # {exc}\n"
                )
                handle.flush()

            processed += 1


def _fragment_atom_list(atoms: str) -> List[int]:
    """Convert an NWChem atom range into a 1-based inclusive atom list."""
    parts = atoms.split(":")
    first = int(parts[0])
    last = int(parts[1])
    return list(range(first, last + 1))


def _fragment_nuclear_charge(atom_types: List[str], fragment_atoms: List[int]) -> float:
    """Return the total nuclear charge (atomic numbers) of a fragment."""
    total = 0.0
    for atom_idx in fragment_atoms:
        symbol = str(atom_types[atom_idx - 1]).capitalize()
        total += _ATOMIC_NUMBERS.get(symbol, 0.0)
    return total
