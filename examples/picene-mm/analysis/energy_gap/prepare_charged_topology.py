#!/usr/bin/env python3
"""Create an Amber topology that differs only in one residue's charges."""

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ChargeChange:
    atom_index_1based: int
    name: str
    neutral_charge: float
    charged_charge: float

    @property
    def delta_charge(self):
        return self.charged_charge - self.neutral_charge


@dataclass(frozen=True)
class TopologyValidationReport:
    status: str
    neutral_target_charge: float
    charged_target_charge: float
    delta_target_charge: float
    target_atom_count: int


def _target_residue(structure, resid_1based):
    if resid_1based < 1 or resid_1based > len(structure.residues):
        raise ValueError(
            f"Residue {resid_1based} is outside the valid 1-based range "
            f"1..{len(structure.residues)}"
        )
    return structure.residues[resid_1based - 1]


def _unique_atom_map(atoms, label):
    mapping = {}
    for item in atoms:
        if item.name in mapping:
            raise ValueError(
                f"Atom names must be unique for name-based mapping; "
                f"duplicate {item.name!r} in {label}"
            )
        mapping[item.name] = item
    return mapping


def transfer_charges_by_name(neutral_structure, charged_mol2, resid_1based):
    """Transfer only atom charges from a MOL2 structure to one residue."""
    target = _target_residue(neutral_structure, resid_1based)
    target_by_name = _unique_atom_map(target.atoms, "target residue")
    mol2_by_name = _unique_atom_map(charged_mol2.atoms, "charged MOL2")

    if len(target.atoms) != len(charged_mol2.atoms):
        raise ValueError(
            "Atom count mismatch: target residue has "
            f"{len(target.atoms)} atoms but charged MOL2 has "
            f"{len(charged_mol2.atoms)}"
        )
    target_names = set(target_by_name)
    mol2_names = set(mol2_by_name)
    if target_names != mol2_names:
        missing = sorted(target_names - mol2_names)
        additional = sorted(mol2_names - target_names)
        raise ValueError(
            "Charged MOL2 atom names do not match the target residue; "
            f"missing={missing}, additional={additional}"
        )

    for name, target_atom in target_by_name.items():
        mol2_atom = mol2_by_name[name]
        target_number = int(getattr(target_atom, "atomic_number", 0) or 0)
        mol2_number = int(getattr(mol2_atom, "atomic_number", 0) or 0)
        if target_number and mol2_number and target_number != mol2_number:
            raise ValueError(
                f"atomic number mismatch for {name}: target={target_number}, "
                f"charged MOL2={mol2_number}"
            )

    changes = []
    for target_atom in target.atoms:
        mol2_atom = mol2_by_name[target_atom.name]
        neutral_charge = float(target_atom.charge)
        charged_charge = float(mol2_atom.charge)
        changes.append(
            ChargeChange(
                atom_index_1based=int(target_atom.idx) + 1,
                name=target_atom.name,
                neutral_charge=neutral_charge,
                charged_charge=charged_charge,
            )
        )
        target_atom.charge = charged_charge
    return changes


def _same_values(left, right):
    try:
        return bool(np.array_equal(np.asarray(left), np.asarray(right)))
    except (TypeError, ValueError):
        return left == right


def validate_topologies(neutral, charged, resid_1based, charge_tolerance=1.0e-8):
    """Fail unless target charges are the only difference between topologies."""
    neutral_target = _target_residue(neutral, resid_1based)
    charged_target = _target_residue(charged, resid_1based)

    if len(neutral.atoms) != len(charged.atoms):
        raise ValueError("Atom counts differ between neutral and charged topologies")
    if len(neutral.residues) != len(charged.residues):
        raise ValueError("Residue counts differ between neutral and charged topologies")

    for index, (left, right) in enumerate(zip(neutral.residues, charged.residues), 1):
        if left.name != right.name or len(left.atoms) != len(right.atoms):
            raise ValueError(f"Residue ordering differs at residue {index}")

    for index, (left, right) in enumerate(zip(neutral.atoms, charged.atoms), 1):
        for attribute in ("name", "type", "mass", "atomic_number"):
            if getattr(left, attribute, None) != getattr(right, attribute, None):
                raise ValueError(
                    f"Atom ordering or parameters differ at atom {index}: {attribute}"
                )

    neutral_flags = set(neutral.parm_data)
    charged_flags = set(charged.parm_data)
    if neutral_flags != charged_flags:
        raise ValueError(
            "Topology data flags differ: "
            f"neutral-only={sorted(neutral_flags - charged_flags)}, "
            f"charged-only={sorted(charged_flags - neutral_flags)}"
        )
    for flag in sorted(neutral_flags - {"CHARGE"}):
        if not _same_values(neutral.parm_data[flag], charged.parm_data[flag]):
            raise ValueError(f"Non-charge topology parameter differs in FLAG {flag}")

    target_indices = {int(item.idx) for item in neutral_target.atoms}
    changed_target_atoms = 0
    for index, (left, right) in enumerate(zip(neutral.atoms, charged.atoms)):
        difference = float(right.charge) - float(left.charge)
        if index not in target_indices and abs(difference) > charge_tolerance:
            raise ValueError(
                f"Non-target charge differs at 1-based atom {index + 1}: "
                f"{left.charge} -> {right.charge}"
            )
        if index in target_indices and abs(difference) > charge_tolerance:
            changed_target_atoms += 1

    if len(neutral_target.atoms) != len(charged_target.atoms):
        raise ValueError("Target residue atom counts differ")
    if changed_target_atoms == 0:
        raise ValueError("Target charges did not change")

    neutral_charge = sum(float(item.charge) for item in neutral_target.atoms)
    charged_charge = sum(float(item.charge) for item in charged_target.atoms)
    return TopologyValidationReport(
        status="PASS",
        neutral_target_charge=neutral_charge,
        charged_target_charge=charged_charge,
        delta_target_charge=charged_charge - neutral_charge,
        target_atom_count=len(neutral_target.atoms),
    )


def print_charge_table(changes):
    print("index   atom       q_neutral       q_charged         delta_q")
    for item in changes:
        print(
            f"{item.atom_index_1based:5d}   {item.name:<6s} "
            f"{item.neutral_charge:15.8f} {item.charged_charge:15.8f} "
            f"{item.delta_charge:15.8f}"
        )


def print_validation_report(report):
    print("\nTopology validation")
    print("-------------------")
    print("Atoms/residues/order:     identical")
    print("Non-charge topology data: identical")
    print("Non-target charges:       identical")
    print("Target charges:           changed")
    print(
        "Target total charge:      "
        f"{report.neutral_target_charge:.8f} -> "
        f"{report.charged_target_charge:.8f} "
        f"(delta {report.delta_target_charge:+.8f})"
    )
    print(f"Status:                   {report.status}")


def load_parmed(path):
    try:
        import parmed
    except ImportError as exc:
        raise RuntimeError(
            "ParmEd is required. Run this script with `conda run -n parmed python ...`."
        ) from exc
    return parmed.load_file(str(path))


def prepare_charged_topology(
    prmtop,
    charged_mol2,
    resid_1based,
    output,
    expected_delta_charge=None,
    charge_tolerance=1.0e-4,
):
    prmtop = Path(prmtop)
    output = Path(output)
    if prmtop.resolve() == output.resolve():
        raise ValueError(
            "Refusing to overwrite the neutral topology; choose a different --output"
        )
    neutral = load_parmed(prmtop)
    charged = load_parmed(prmtop)
    mol2 = load_parmed(charged_mol2)
    changes = transfer_charges_by_name(charged, mol2, resid_1based)
    print_charge_table(changes)

    output.parent.mkdir(parents=True, exist_ok=True)
    charged.save(str(output), format="amber", overwrite=True)
    written = load_parmed(output)
    report = validate_topologies(neutral, written, resid_1based)
    print_validation_report(report)

    if expected_delta_charge is not None and not np.isclose(
        report.delta_target_charge,
        expected_delta_charge,
        atol=charge_tolerance,
        rtol=0.0,
    ):
        raise ValueError(
            f"Target charge change {report.delta_target_charge:+.8f} does not match "
            f"expected {expected_delta_charge:+.8f} within {charge_tolerance:g}"
        )
    return report


def build_parser():
    parser = argparse.ArgumentParser(
        description="Replace only one Amber residue's charges using a charged-state MOL2"
    )
    parser.add_argument("--prmtop", required=True, type=Path)
    parser.add_argument("--charged-mol2", required=True, type=Path)
    parser.add_argument(
        "--resid", required=True, type=int, help="Target residue number (1-based)"
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-delta-charge", type=float)
    parser.add_argument("--charge-tolerance", type=float, default=1.0e-4)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    prepare_charged_topology(
        args.prmtop,
        args.charged_mol2,
        args.resid,
        args.output,
        expected_delta_charge=args.expected_delta_charge,
        charge_tolerance=args.charge_tolerance,
    )


if __name__ == "__main__":
    main()
