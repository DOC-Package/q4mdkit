"""Preserve selected OpenMM constraints during QM/MM theory construction."""


def capture_qm_hbond_constraints(system, qmatoms, elements):
    """Return constraints on covalent X-H bonds involving a QM atom.

    Here HBonds has OpenMM's meaning: bonds to hydrogen, rather than
    intermolecular hydrogen bonds. Distances retain their original units.
    """
    qm_atoms = set(qmatoms)
    constraints = []
    for index in range(system.getNumConstraints()):
        atom1, atom2, distance = system.getConstraintParameters(index)
        atom1, atom2 = int(atom1), int(atom2)
        if ((atom1 in qm_atoms or atom2 in qm_atoms)
                and (elements[atom1].upper() == 'H' or elements[atom2].upper() == 'H')):
            constraints.append((atom1, atom2, distance))
    return constraints


def restore_constraints(system, constraints):
    """Restore captured constraints and return how many were added."""
    for atom1, atom2, distance in constraints:
        system.addConstraint(atom1, atom2, distance)
    return len(constraints)
