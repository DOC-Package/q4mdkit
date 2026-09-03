#!/usr/bin/env python3
"""Evaluate a vertical neutral-to-charged MM energy gap on a trajectory."""

import argparse
import csv
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


AVOGADRO_CONSTANT_MOL = 6.02214076e23
PLANCK_CONSTANT_J_S = 6.62607015e-34
SPEED_OF_LIGHT_M_S = 299792458.0
KJMOL_TO_CM1 = 1000.0 / (
    AVOGADRO_CONSTANT_MOL * PLANCK_CONSTANT_J_S * SPEED_OF_LIGHT_M_S * 100.0
)
AMBER_CHARGE_SCALE = 18.2223


@dataclass(frozen=True)
class EnergyEvaluation:
    frame: int
    time_ps: float
    neutral_kjmol: float
    charged_kjmol: float

    @property
    def gap_kjmol(self):
        return self.charged_kjmol - self.neutral_kjmol


@dataclass(frozen=True)
class AmberTopologyValidation:
    changed_target_atoms: int
    delta_target_charge: float


def validate_amber_topology_data(
    neutral_loader, charged_loader, resid_1based, charge_tolerance=1.0e-8
):
    """Check OpenMM-loaded prmtop data before evaluating an energy gap."""
    neutral_data = neutral_loader._raw_data
    charged_data = charged_loader._raw_data
    neutral_flags = set(neutral_data)
    charged_flags = set(charged_data)
    if neutral_flags != charged_flags:
        raise ValueError(
            "Amber topology flags differ: "
            f"neutral-only={sorted(neutral_flags - charged_flags)}, "
            f"charged-only={sorted(charged_flags - neutral_flags)}"
        )
    for flag in sorted(neutral_flags - {"CHARGE"}):
        if neutral_data[flag] != charged_data[flag]:
            raise ValueError(f"Non-charge Amber topology FLAG differs: {flag}")

    residue_pointers = [int(value) - 1 for value in neutral_data["RESIDUE_POINTER"]]
    if resid_1based < 1 or resid_1based > len(residue_pointers):
        raise ValueError(
            f"Residue {resid_1based} is outside the valid 1-based range "
            f"1..{len(residue_pointers)}"
        )
    atom_count = len(neutral_data["ATOM_NAME"])
    target_start = residue_pointers[resid_1based - 1]
    target_stop = (
        residue_pointers[resid_1based]
        if resid_1based < len(residue_pointers)
        else atom_count
    )
    neutral_charges = np.asarray(neutral_data["CHARGE"], dtype=float) / AMBER_CHARGE_SCALE
    charged_charges = np.asarray(charged_data["CHARGE"], dtype=float) / AMBER_CHARGE_SCALE
    if neutral_charges.shape != charged_charges.shape or len(neutral_charges) != atom_count:
        raise ValueError("Amber topology charge arrays are incompatible")

    changed_target_atoms = 0
    for index, difference in enumerate(charged_charges - neutral_charges):
        in_target = target_start <= index < target_stop
        if not in_target and abs(float(difference)) > charge_tolerance:
            raise ValueError(
                f"Non-target charge differs at 1-based atom {index + 1}: "
                f"{neutral_charges[index]} -> {charged_charges[index]}"
            )
        if in_target and abs(float(difference)) > charge_tolerance:
            changed_target_atoms += 1
    delta_target_charge = float(
        np.sum(charged_charges[target_start:target_stop])
        - np.sum(neutral_charges[target_start:target_stop])
    )
    return AmberTopologyValidation(
        changed_target_atoms=changed_target_atoms,
        delta_target_charge=delta_target_charge,
    )


def center_energy_gaps(gaps):
    values = np.asarray(gaps, dtype=float)
    if values.size == 0:
        raise ValueError("Cannot center an empty energy-gap trajectory")
    mean = float(np.mean(values))
    return mean, (values - mean).tolist()


def frame_time_ps(frame_index, reader_time_ps, origin_ps=None, step_ps=None):
    if (origin_ps is None) != (step_ps is None):
        raise ValueError("--time-origin-ps and --time-step-ps must be supplied together")
    if origin_ps is not None:
        return float(origin_ps + frame_index * step_ps)
    return float(reader_time_ps)


def write_results_csv(output, evaluations):
    evaluations = list(evaluations)
    mean_gap, centered = center_energy_gaps(
        [item.gap_kjmol for item in evaluations]
    )
    fieldnames = [
        "frame",
        "time_ps",
        "E_neutral_kJmol",
        "E_charged_kJmol",
        "gap_kJmol",
        "gap_fluctuation_kJmol",
        "gap_cm-1",
        "gap_fluctuation_cm-1",
    ]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item, fluctuation in zip(evaluations, centered):
            gap = item.gap_kjmol
            writer.writerow(
                {
                    "frame": item.frame,
                    "time_ps": f"{item.time_ps:.12g}",
                    "E_neutral_kJmol": f"{item.neutral_kjmol:.12g}",
                    "E_charged_kJmol": f"{item.charged_kjmol:.12g}",
                    "gap_kJmol": f"{gap:.12g}",
                    "gap_fluctuation_kJmol": f"{fluctuation:.12g}",
                    "gap_cm-1": f"{gap * KJMOL_TO_CM1:.12g}",
                    "gap_fluctuation_cm-1": f"{fluctuation * KJMOL_TO_CM1:.12g}",
                }
            )
    return mean_gap, centered


def _openmm_modules():
    try:
        import openmm
        from openmm import app, unit
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM is required. Run this script with `conda run -n myash python ...`."
        ) from exc
    return openmm, app, unit


def configure_periodic_amber_topology(
    amber_prmtop, box_vectors_nm, nanometer_unit
):
    """Attach a box and enable PBC for Amber files written with IFBOX=0."""
    vectors = np.asarray(box_vectors_nm, dtype=float) * nanometer_unit
    amber_prmtop.topology.setPeriodicBoxVectors(vectors)
    try:
        amber_prmtop._prmtop._raw_data["POINTERS"][27] = 1
    except (AttributeError, KeyError, IndexError) as exc:
        raise ValueError("Could not enable IFBOX in the Amber topology") from exc


def load_initial_box_vectors_nm(trajectory, topology, start):
    try:
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError(
            "MDTraj is required. Run this script with `conda run -n myash python ...`."
        ) from exc
    try:
        frame = md.load_frame(str(trajectory), index=start, top=str(topology))
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Could not read trajectory frame {start}") from exc
    if frame.unitcell_vectors is None:
        raise ValueError(f"Periodic box vectors are missing at frame {start}")
    return np.asarray(frame.unitcell_vectors[0], dtype=float)


def create_openmm_system(
    prmtop_path,
    nonbonded_method="PME",
    cutoff_angstrom=9.0,
    ewald_error_tolerance=5.0e-4,
    initial_box_vectors_nm=None,
):
    openmm, app, unit = _openmm_modules()
    del openmm
    method_map = {
        "PME": app.PME,
        "Ewald": app.Ewald,
        "CutoffPeriodic": app.CutoffPeriodic,
        "NoCutoff": app.NoCutoff,
    }
    topology = app.AmberPrmtopFile(str(prmtop_path))
    if nonbonded_method != "NoCutoff":
        if initial_box_vectors_nm is None:
            raise ValueError(
                "Periodic system construction requires initial box vectors"
            )
        configure_periodic_amber_topology(
            topology, initial_box_vectors_nm, unit.nanometer
        )
    kwargs = {
        "nonbondedMethod": method_map[nonbonded_method],
        "constraints": app.HBonds,
        "rigidWater": False,
        "hydrogenMass": 1.0 * unit.dalton,
    }
    if nonbonded_method != "NoCutoff":
        kwargs.update(
            nonbondedCutoff=cutoff_angstrom * unit.angstrom,
            ewaldErrorTolerance=ewald_error_tolerance,
        )
    system = topology.createSystem(**kwargs)
    for force in system.getForces():
        if force.__class__.__name__ == "NonbondedForce":
            force.setUseDispersionCorrection(True)
    return system


def validate_amber_topology_files(neutral_prmtop, charged_prmtop, resid_1based):
    openmm, app, unit = _openmm_modules()
    del openmm, unit
    neutral = app.AmberPrmtopFile(str(neutral_prmtop))
    charged = app.AmberPrmtopFile(str(charged_prmtop))
    return validate_amber_topology_data(
        neutral._prmtop, charged._prmtop, resid_1based
    )


def _platform_and_properties(platform_name, precision, threads):
    openmm, app, unit = _openmm_modules()
    del app, unit
    platform = openmm.Platform.getPlatformByName(platform_name)
    supported = set(platform.getPropertyNames())
    properties = {}
    if "Threads" in supported:
        properties["Threads"] = str(threads)
    if precision != "default" and "Precision" in supported:
        properties["Precision"] = precision
    return platform, properties


def _create_context(system, platform, properties):
    openmm, app, unit = _openmm_modules()
    del app
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    context = openmm.Context(system, integrator, platform, properties)
    return context, integrator


def _selected_frames(trajectory, topology, start, stop, stride, chunk_size):
    try:
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError(
            "MDTraj is required. Run this script with `conda run -n myash python ...`."
        ) from exc
    if start < 0:
        raise ValueError("--start must be non-negative")
    if stride < 1:
        raise ValueError("--stride must be at least 1")
    if stop is not None and stop <= start:
        raise ValueError("--stop must be greater than --start")
    selected_index = 0
    for chunk in md.iterload(
        str(trajectory), top=str(topology), chunk=chunk_size, skip=start, stride=stride
    ):
        for local_index in range(chunk.n_frames):
            frame_index = start + selected_index * stride
            selected_index += 1
            if stop is not None and frame_index >= stop:
                return
            yield frame_index, chunk[local_index]


def evaluate_energy_gap(
    neutral_prmtop,
    charged_prmtop,
    trajectory,
    resid_1based,
    start=0,
    stop=None,
    stride=1,
    time_origin_ps=None,
    time_step_ps=None,
    platform_name="CPU",
    precision="default",
    threads=4,
    nonbonded_method="PME",
    cutoff_angstrom=9.0,
    ewald_error_tolerance=5.0e-4,
    chunk_size=100,
):
    openmm, app, unit = _openmm_modules()
    del app
    topology_report = validate_amber_topology_files(
        neutral_prmtop, charged_prmtop, resid_1based
    )
    print(
        "Topology guard: non-charge FLAGs and non-target charges identical; "
        f"changed target atoms={topology_report.changed_target_atoms}, "
        f"target delta charge={topology_report.delta_target_charge:+.8f}"
    )
    initial_box_vectors_nm = None
    if nonbonded_method != "NoCutoff":
        initial_box_vectors_nm = load_initial_box_vectors_nm(
            trajectory, neutral_prmtop, start
        )
    neutral_system = create_openmm_system(
        neutral_prmtop,
        nonbonded_method,
        cutoff_angstrom,
        ewald_error_tolerance,
        initial_box_vectors_nm,
    )
    charged_system = create_openmm_system(
        charged_prmtop,
        nonbonded_method,
        cutoff_angstrom,
        ewald_error_tolerance,
        initial_box_vectors_nm,
    )
    print("Created neutral and charged OpenMM systems.", flush=True)
    if neutral_system.getNumParticles() != charged_system.getNumParticles():
        raise ValueError("Neutral and charged OpenMM systems have different atom counts")

    platform, properties = _platform_and_properties(platform_name, precision, threads)
    neutral_context, neutral_integrator = _create_context(
        neutral_system, platform, properties
    )
    charged_context, charged_integrator = _create_context(
        charged_system, platform, properties
    )
    del neutral_integrator, charged_integrator
    print(
        f"Created OpenMM contexts on {platform_name}; evaluating selected frames...",
        flush=True,
    )

    evaluations = []
    total_frames = None
    if stop is not None:
        total_frames = max(0, math.ceil((stop - start) / stride))
    processed = 0
    start_time = time.perf_counter()
    for frame_index, frame in _selected_frames(
        trajectory, neutral_prmtop, start, stop, stride, chunk_size
    ):
        if frame.n_atoms != neutral_system.getNumParticles():
            raise ValueError(
                f"Trajectory frame has {frame.n_atoms} atoms but topology has "
                f"{neutral_system.getNumParticles()}"
            )
        if nonbonded_method != "NoCutoff":
            if frame.unitcell_vectors is None:
                raise ValueError(f"Periodic box vectors are missing at frame {frame_index}")
            vectors = [
                openmm.Vec3(*map(float, vector)) * unit.nanometer
                for vector in frame.unitcell_vectors[0]
            ]
            neutral_context.setPeriodicBoxVectors(*vectors)
            charged_context.setPeriodicBoxVectors(*vectors)
        positions = frame.xyz[0] * unit.nanometer
        neutral_context.setPositions(positions)
        charged_context.setPositions(positions)
        neutral_energy = neutral_context.getState(getEnergy=True).getPotentialEnergy()
        charged_energy = charged_context.getState(getEnergy=True).getPotentialEnergy()
        neutral_value = float(neutral_energy.value_in_unit(unit.kilojoule_per_mole))
        charged_value = float(charged_energy.value_in_unit(unit.kilojoule_per_mole))
        if not math.isfinite(neutral_value) or not math.isfinite(charged_value):
            raise ValueError(f"Non-finite energy at frame {frame_index}")
        evaluations.append(
            EnergyEvaluation(
                frame=frame_index,
                time_ps=frame_time_ps(
                    frame_index,
                    float(frame.time[0]),
                    origin_ps=time_origin_ps,
                    step_ps=time_step_ps,
                ),
                neutral_kjmol=neutral_value,
                charged_kjmol=charged_value,
            )
        )
        processed += 1
        elapsed = time.perf_counter() - start_time
        progress_count = f"{processed}/{total_frames}" if total_frames else str(processed)
        print(
            f"Progress [{progress_count}] frame={frame_index} "
            f"time_ps={evaluations[-1].time_ps:.6f} "
            f"neutral={neutral_value:.8f} charged={charged_value:.8f} "
            f"gap={charged_value - neutral_value:+.8f} kJ/mol "
            f"elapsed={elapsed:.1f}s",
            flush=True,
        )
    if not evaluations:
        raise ValueError("Frame selection produced no trajectory frames")
    print(
        f"Evaluated {len(evaluations)} frames in "
        f"{time.perf_counter() - start_time:.1f}s",
        flush=True,
    )
    return evaluations, properties


def print_summary(evaluations, platform_name, precision, properties):
    gaps = np.asarray([item.gap_kjmol for item in evaluations], dtype=float)
    print("\nEnergy-gap summary")
    print("------------------")
    print(f"OpenMM platform:   {platform_name}")
    print(f"Requested precision: {precision}")
    print(f"Applied properties: {properties or 'platform defaults'}")
    print(f"Number of frames:  {len(gaps)}")
    print(f"Mean energy gap:   {np.mean(gaps):.10f} kJ/mol")
    print(f"Standard deviation:{np.std(gaps):.10f} kJ/mol")
    print(f"Minimum:           {np.min(gaps):.10f} kJ/mol")
    print(f"Maximum:           {np.max(gaps):.10f} kJ/mol")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Compute U_charged(R_t)-U_neutral(R_t) on an existing trajectory"
    )
    parser.add_argument("--neutral-prmtop", required=True, type=Path)
    parser.add_argument("--charged-prmtop", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--resid", required=True, type=int, help="Target residue number (1-based)"
    )
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--time-origin-ps", type=float)
    parser.add_argument("--time-step-ps", type=float)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument(
        "--precision", choices=("default", "single", "mixed", "double"), default="default"
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--nonbonded-method",
        choices=("PME", "Ewald", "CutoffPeriodic", "NoCutoff"),
        default="PME",
    )
    parser.add_argument("--cutoff-angstrom", type=float, default=9.0)
    parser.add_argument("--ewald-error-tolerance", type=float, default=5.0e-4)
    parser.add_argument("--chunk-size", type=int, default=100)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be at least 1")
    print("Vertical energy-gap evaluation (no minimization or dynamics)")
    print(f"Neutral topology:  {args.neutral_prmtop}")
    print(f"Charged topology:  {args.charged_prmtop}")
    print(f"Trajectory:        {args.trajectory}")
    print(
        f"OpenMM settings:   method={args.nonbonded_method}, "
        f"cutoff={args.cutoff_angstrom:g} A, "
        f"Ewald tolerance={args.ewald_error_tolerance:g}"
    )
    evaluations, properties = evaluate_energy_gap(
        args.neutral_prmtop,
        args.charged_prmtop,
        args.trajectory,
        args.resid,
        start=args.start,
        stop=args.stop,
        stride=args.stride,
        time_origin_ps=args.time_origin_ps,
        time_step_ps=args.time_step_ps,
        platform_name=args.platform,
        precision=args.precision,
        threads=args.threads,
        nonbonded_method=args.nonbonded_method,
        cutoff_angstrom=args.cutoff_angstrom,
        ewald_error_tolerance=args.ewald_error_tolerance,
        chunk_size=args.chunk_size,
    )
    write_results_csv(args.output, evaluations)
    print_summary(evaluations, args.platform, args.precision, properties)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
