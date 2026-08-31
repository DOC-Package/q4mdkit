"""Helpers for preserving selected OpenMM constraints across QMMM setup."""


def capture_qm_hbond_constraints(system, qmatoms, elements):
    """Return X-H constraints involving at least one QM atom."""
    qm_atom_set = set(qmatoms)
    constraints = []
    for index in range(system.getNumConstraints()):
        atom1, atom2, distance = system.getConstraintParameters(index)
        atom1 = int(atom1)
        atom2 = int(atom2)
        involves_qm = atom1 in qm_atom_set or atom2 in qm_atom_set
        involves_hydrogen = elements[atom1].upper() == "H" or elements[atom2].upper() == "H"
        if involves_qm and involves_hydrogen:
            constraints.append((atom1, atom2, distance))
    return constraints


def restore_constraints(system, constraints):
    """Add previously captured constraints to an OpenMM System."""
    for atom1, atom2, distance in constraints:
        system.addConstraint(atom1, atom2, distance)
    return len(constraints)
