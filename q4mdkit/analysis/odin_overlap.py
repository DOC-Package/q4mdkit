"""
Interface to the ODIN program (T. Niehaus) for computing exact cross-geometry
AO overlap matrices between consecutive MD frames.

ODIN computes the DFTB overlap matrix for an arbitrary geometry using
Slater-Koster tables and writes the result to ``oversqr.dat``.

Cross-geometry strategy
-----------------------
The exact cross-geometry AO overlap between frames t_{n-1} and t_n is

    S^{n-1,n}_{mu,nu} = <chi_mu(t_{n-1}) | chi_nu(t_n)>

where basis functions are centred on atoms at *different* geometries.
We compute this by the *doubled-geometry trick*:

1.  Create a gen file with 2N atoms: atoms 1..N at positions r(t_{n-1}),
    atoms N+1..2N at positions r(t_n), all with their real element labels.
2.  Run ODIN on this 2N-atom system.  Because DFTB Slater-Koster integrals
    depend only on the inter-atomic displacement vector and element types,
    the off-diagonal block of the (2*N_AO x 2*N_AO) output matrix gives
    exactly the desired cross-geometry overlaps.
3.  Extract the upper-left / lower-right cross block of size N_AO x N_AO.

Note on special SK files
-------------------------
For short inter-atomic distances (which can appear in cross-geometry pairs
even when both individual geometries are well-behaved), the standard SK files
from dftb.org contain only dummy values (typically 20 * 1.0).  ODIN handles
this gracefully, but accurate results require SK files generated with
GridStart = GridSeparation (e.g., via the skprogs toolchain).  See the ODIN
README for details.

Note on AO ordering
--------------------
ODIN internally reorders AOs via the ``icv`` index array.  The ordering used
here matches the DFTB+ convention for real-space cluster calculations, so the
cross-overlap matrix can be used directly with MO coefficient matrices
produced by DFTB+.  If you suspect a mismatch, compare the diagonal block of
the doubled-geometry overlap with the single-geometry S_AO from DFTB+.
"""

from __future__ import annotations

import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


_ATOMIC_MASSES = {
    'H': 1.008,
    'He': 4.002602,
    'Li': 6.941,
    'Be': 9.012182,
    'B': 10.811,
    'C': 12.01,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984032,
    'Ne': 20.1797,
    'Na': 22.98976928,
    'Mg': 24.3050,
    'Al': 26.9815386,
    'Si': 28.0855,
    'P': 30.973762,
    'S': 32.065,
    'Cl': 35.453,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cross_overlap_odin(
    coords_prev_ang: np.ndarray,
    coords_curr_ang: np.ndarray,
    atom_types_per_atom: List[str],
    lmax_dict: Dict[str, int],
    sk_prefix: str,
    sk_separator: str,
    sk_suffix: str,
    odin_executable: str,
    remove_translation_rotation: bool = False,
    work_dir: Optional[Path] = None,
    keep_files: bool = False,
) -> np.ndarray:
    """
    Compute the cross-geometry AO overlap matrix S^{n-1,n} using ODIN.

    Parameters
    ----------
    coords_prev_ang : np.ndarray, shape (N, 3)
        QM-region Cartesian coordinates at t_{n-1}, in Angstrom.
    coords_curr_ang : np.ndarray, shape (N, 3)
        QM-region Cartesian coordinates at t_n, in Angstrom.
    atom_types_per_atom : list of str
        Element symbol for each of the N atoms (same order as ``coords_*``).
    lmax_dict : dict {element: int}
        Maximum angular-momentum quantum number per element in the ODIN
        convention: 1=s, 2=p, 3=d.  Must cover all elements in
        ``atom_types_per_atom``.
    sk_prefix : str
        Path prefix for SK files, e.g. ``'/home/user/sk/3ob/'``.
    sk_separator : str
        Separator between element symbols in SK file names, e.g. ``'-'``.
    sk_suffix : str
        Suffix for SK files, e.g. ``'.skf'``.
    odin_executable : str
        Absolute path to the compiled ``odin`` binary.
    remove_translation_rotation : bool
        If True, remove rigid translation and rotation from the two geometries
        before building the doubled ODIN system. This follows the same
        center-of-mass removal plus mass-weighted rigid-body alignment used in
        ``normal_mode_analysis.py``.
    work_dir : Path, optional
        Directory in which to run ODIN.  A temporary directory is used when
        *None* and cleaned up afterwards (unless ``keep_files=True``).
    keep_files : bool
        If True, preserve all ODIN inputs/outputs for each frame in the
        caller-provided work directory. If False and ``work_dir`` is provided,
        the latest ODIN files are kept in that directory and overwritten on the
        next frame. Temporary auto-created directories are still removed.

    Returns
    -------
    S_cross : np.ndarray, shape (N_AO, N_AO)
        Cross-geometry AO overlap matrix where row index μ refers to AOs at
        t_{n-1} and column index ν refers to AOs at t_n.
    """
    N = len(atom_types_per_atom)
    if coords_prev_ang.shape != (N, 3) or coords_curr_ang.shape != (N, 3):
        raise ValueError("coords_prev_ang and coords_curr_ang must both be (N, 3)")

    if remove_translation_rotation:
        coords_prev_ang, coords_curr_ang = _remove_translation_rotation_from_geometries(
            coords_prev_ang,
            coords_curr_ang,
            atom_types_per_atom,
        )

    # Unique element types in gen-file order (order of first appearance)
    unique_types = _unique_ordered(atom_types_per_atom)
    lmax_list = [lmax_dict[t] for t in unique_types]

    # Number of AOs per atom and total for one geometry
    n_ao_per_atom = [lmax_dict[t] ** 2 for t in atom_types_per_atom]
    N_AO = sum(n_ao_per_atom)

    # Set up working directory
    tmp_dir: Optional[str] = None
    if work_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix="odin_cross_")
        run_dir = Path(tmp_dir)
    else:
        run_dir = Path(work_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

    try:
        gen_path = run_dir / "cross_geo.gen"
        inp_path = run_dir / "odin.inp"
        out_path = run_dir / "oversqr.dat"

        _write_doubled_gen(
            coords_prev_ang, coords_curr_ang, atom_types_per_atom, gen_path
        )
        _write_odin_inp(gen_path, sk_prefix, sk_separator, sk_suffix, lmax_list, inp_path)
        _run_odin(odin_executable, inp_path, run_dir)

        full_matrix = _read_oversqr(out_path)
        if full_matrix.shape[0] != 2 * N_AO:
            raise RuntimeError(
                f"ODIN returned matrix of size {full_matrix.shape[0]}, "
                f"expected {2 * N_AO} (= 2 x N_AO for doubled system)."
            )

        # Upper-right block: rows = t_{n-1} AOs, cols = t_n AOs
        S_cross = full_matrix[:N_AO, N_AO:].copy()
        return S_cross

    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def n_ao_for_atoms(atom_types_per_atom: List[str], lmax_dict: Dict[str, int]) -> int:
    """Return the total number of AOs for a list of atoms given lmax values."""
    return sum(lmax_dict[t] ** 2 for t in atom_types_per_atom)


def _remove_translation_rotation_from_geometries(
    coords_prev_ang: np.ndarray,
    coords_curr_ang: np.ndarray,
    atom_types_per_atom: List[str],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Remove rigid translation and rotation before cross-overlap evaluation.

    The previous geometry is shifted to its center-of-mass frame. The current
    geometry is shifted to its own center of mass and then mass-weight aligned
    onto the previous geometry using the same rigid-body-removal procedure as
    in ``normal_mode_analysis.py``.
    """
    masses = np.array([_ATOMIC_MASSES[t] for t in atom_types_per_atom], dtype=float)
    prev_centered = _remove_center_of_mass(coords_prev_ang, masses)
    curr_centered = _remove_center_of_mass(coords_curr_ang, masses)
    rotation = _weighted_kabsch_rotation(curr_centered, prev_centered, masses)
    curr_aligned = curr_centered @ rotation
    return prev_centered, curr_aligned


def _remove_center_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Return coordinates centered at the mass-weighted center of mass."""
    total_mass = float(np.sum(masses))
    center = np.sum(masses[:, np.newaxis] * coords, axis=0) / total_mass
    return coords - center


def _weighted_kabsch_rotation(
    source: np.ndarray,
    target: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Return the mass-weighted rigid-body rotation aligning ``source`` onto ``target``."""
    covariance = (weights[:, np.newaxis] * source).T @ target
    U, _, Vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(U @ Vt) < 0.0:
        correction[-1, -1] = -1.0
    return U @ correction @ Vt


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _unique_ordered(items: List[str]) -> List[str]:
    """Return unique items in the order of their first appearance."""
    seen: List[str] = []
    for x in items:
        if x not in seen:
            seen.append(x)
    return seen


def _write_doubled_gen(
    coords_prev_ang: np.ndarray,
    coords_curr_ang: np.ndarray,
    atom_types_per_atom: List[str],
    filepath: Path,
) -> None:
    """
    Write a DFTB+ gen file with 2N atoms:
        atoms 1..N   at coords_prev_ang  (geometry t_{n-1})
        atoms N+1..2N at coords_curr_ang (geometry t_n)
    All atoms keep their real element labels so ODIN uses the correct SK files.
    """
    N = len(atom_types_per_atom)
    unique_types = _unique_ordered(atom_types_per_atom)
    type_idx = {t: i + 1 for i, t in enumerate(unique_types)}

    with open(filepath, "w") as f:
        f.write(f" {2 * N}  C\n")
        f.write("  " + "  ".join(unique_types) + "\n")
        for i, (r, t) in enumerate(zip(coords_prev_ang, atom_types_per_atom)):
            f.write(
                f"  {i + 1}  {type_idx[t]}"
                f"  {r[0]:16.8f}  {r[1]:16.8f}  {r[2]:16.8f}\n"
            )
        for i, (r, t) in enumerate(zip(coords_curr_ang, atom_types_per_atom)):
            f.write(
                f"  {N + i + 1}  {type_idx[t]}"
                f"  {r[0]:16.8f}  {r[1]:16.8f}  {r[2]:16.8f}\n"
            )


def _write_odin_inp(
    gen_path: Path,
    sk_prefix: str,
    sk_sep: str,
    sk_suffix: str,
    lmax_list: List[int],
    inp_path: Path,
) -> None:
    """
    Write the ODIN stdin input file.

    ODIN reads five items interactively:
        1. geometry filename (gen format)
        2. SK file prefix
        3. SK file separator (between element symbols)
        4. SK file suffix
        5. lmax values (integers) for each element type in gen-file order
    """
    lmax_str = "  ".join(str(l) for l in lmax_list)
    with open(inp_path, "w") as f:
        # ODIN uses Fortran list-directed input (read *,...).
        # Single-quoted strings are safe across gfortran versions.
        f.write(f"'{gen_path}'\n")
        f.write(f"'{sk_prefix}'\n")
        f.write(f"'{sk_sep}'\n")
        f.write(f"'{sk_suffix}'\n")
        f.write(f"{lmax_str}\n")


def _run_odin(executable: str, inp_path: Path, run_dir: Path) -> None:
    """Run ODIN with stdin redirected from ``inp_path``."""
    with open(inp_path) as inp_fh:
        result = subprocess.run(
            [executable],
            stdin=inp_fh,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(run_dir),
        )
    if result.returncode != 0:
        stderr_msg = result.stderr.decode(errors="replace")
        stdout_msg = result.stdout.decode(errors="replace")
        raise RuntimeError(
            f"ODIN exited with code {result.returncode}.\n"
            f"stdout:\n{stdout_msg}\nstderr:\n{stderr_msg}"
        )


def _read_oversqr(filepath: Path) -> np.ndarray:
    """
    Read ODIN's ``oversqr.dat`` and return the full overlap matrix.

    File format (written by ODIN in Fortran column-major order)::

        # REAL      NALLORB  NKPOINT
          T          ndim     1
        # IKPOINT
          1
        # MATRIX
          val1  val2  ...  (ndim values per line, column-major order)

    The matrix elements satisfy ``overl(row, col)`` in Fortran (1-indexed),
    so after reading with Fortran reshape order the Python result is
    ``result[row-1, col-1]``.
    """
    with open(filepath) as fh:
        lines = fh.readlines()

    # Line index 1 (0-based): "  T      ndim    1"
    ndim = int(lines[1].split()[1])

    # Matrix data starts at line 5 (after two header lines, ikpoint header,
    # ikpoint value, matrix header).
    values: List[float] = []
    for line in lines[5:]:
        values.extend(float(x) for x in line.split())

    if len(values) != ndim * ndim:
        raise ValueError(
            f"oversqr.dat: expected {ndim * ndim} matrix elements, "
            f"got {len(values)}."
        )

    # Reshape with Fortran (column-major) order so that result[i,j]
    # corresponds to overl(i+1, j+1) in Fortran 1-indexed convention.
    return np.array(values, dtype=float).reshape((ndim, ndim), order="F")
