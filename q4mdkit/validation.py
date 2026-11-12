"""Dependency-light validation of AMBER MD inputs for arbitrary systems."""

import math
from numbers import Integral
from pathlib import Path


def load_box_vectors(box_path):
    """Read 3x3 vectors or legacy ``a b c alpha beta gamma`` cell parameters.

    Lengths are in angstrom and angles in degrees.
    """
    rows = []
    for line in Path(box_path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        values = [float(value) for value in line.split()]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"invalid box-vector row in {box_path}: {line}")
        rows.append(values)
    if len(rows) == 1 and len(rows[0]) == 6:
        a, b, c, alpha, beta, gamma = rows[0]
        if min(a, b, c) <= 0 or any(not 0 < angle < 180 for angle in (alpha, beta, gamma)):
            raise ValueError(f"invalid cell parameters in {box_path}")
        alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
        cx = c * math.cos(beta)
        cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        cz_squared = c * c - cx * cx - cy * cy
        if cz_squared <= 0:
            raise ValueError(f"degenerate cell parameters in {box_path}")
        return [[a, 0.0, 0.0], [b * math.cos(gamma), b * math.sin(gamma), 0.0],
                [cx, cy, math.sqrt(cz_squared)]]
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        raise ValueError(f"expected three box vectors or six cell parameters in {box_path}")
    return rows


def read_prmtop_atom_count(prmtop_path):
    """Read NATOM from an AMBER topology without importing ASH/OpenMM."""
    lines = Path(prmtop_path).read_text().splitlines()
    try:
        flag = next(index for index, line in enumerate(lines)
                    if line.strip() == '%FLAG POINTERS')
        return int(lines[flag + 2].split()[0])
    except (StopIteration, IndexError, ValueError) as error:
        raise ValueError(f"cannot read atom count from {prmtop_path}") from error


def validate_amber_system(pdb_path, prmtop_path, box_path=None, expected_atoms=None,
                          residue_blocks=None, qm_atoms=None, expected_qm_atoms=None,
                          qm_resname=None, periodic=True):
    """Check MD input consistency and return the number of atoms.

    ``expected_atoms`` and ``residue_blocks`` are optional system-specific
    checks. A residue block is ``{'name': 'SOL', 'count': 100}`` and blocks
    describe consecutive PDB atoms. QM selections may be noncontiguous.
    Box-file lengths and PDB CRYST1 lengths must agree within 0.001 angstrom.
    ``periodic=False`` permits PDB inputs without a CRYST1 record.
    """
    lines = Path(pdb_path).read_text().splitlines()
    atoms = [line for line in lines if line.startswith(('ATOM  ', 'HETATM'))]
    n_atoms = len(atoms)
    if n_atoms == 0:
        raise ValueError(f"PDB has no atoms: {pdb_path}")
    topology_atoms = read_prmtop_atom_count(prmtop_path)
    if topology_atoms != n_atoms:
        raise ValueError(f"topology has {topology_atoms} atoms but PDB has {n_atoms}")
    if expected_atoms is not None and n_atoms != expected_atoms:
        raise ValueError(f"PDB has {n_atoms} atoms; expected {expected_atoms}")

    residues = [line[17:20].strip() for line in atoms]
    if residue_blocks is not None:
        offset = 0
        for block in residue_blocks:
            count = block['count']
            if not isinstance(count, Integral) or count <= 0:
                raise ValueError('residue block counts must be positive integers')
            end = offset + count
            if end > n_atoms or any(name != block['name'] for name in residues[offset:end]):
                raise ValueError(f"PDB residue ordering does not match {block['name']} block at atom {offset}")
            offset = end
        if offset != n_atoms:
            raise ValueError(f"residue blocks describe {offset} atoms but PDB has {n_atoms}")

    if qm_atoms is not None:
        selection = list(qm_atoms)
        if (not selection or any(not isinstance(index, Integral) for index in selection)
                or len(set(selection)) != len(selection)
                or any(index < 0 or index >= n_atoms for index in selection)):
            raise ValueError('QM atom indices must be nonempty, unique integers within the PDB atom range')
        if expected_qm_atoms is not None and len(selection) != expected_qm_atoms:
            raise ValueError(f"QM region has {len(selection)} atoms; expected {expected_qm_atoms}")
        if qm_resname is not None and any(residues[index] != qm_resname for index in selection):
            raise ValueError(f"QM atoms must have residue name {qm_resname}")

    cryst1 = next((line for line in lines if line.startswith('CRYST1')), None)
    if periodic and cryst1 is None:
        raise ValueError(f"PDB has no CRYST1 record: {pdb_path}")
    if box_path is not None:
        if cryst1 is None:
            raise ValueError(f"PDB has no CRYST1 record: {pdb_path}")
        vectors = load_box_vectors(box_path)
        lengths = [math.sqrt(sum(value * value for value in vector)) for vector in vectors]
        pdb_lengths = [float(cryst1[start:end]) for start, end in ((6, 15), (15, 24), (24, 33))]
        if (any(not math.isfinite(value) or value <= 0 for value in lengths + pdb_lengths)
                or any(abs(pdb - box) > 1e-3 for pdb, box in zip(pdb_lengths, lengths))):
            raise ValueError(f"PDB CRYST1 box {pdb_lengths} does not match box file {lengths}")
    return n_atoms
