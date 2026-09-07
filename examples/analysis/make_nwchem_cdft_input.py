#!/usr/bin/env python
"""Generate a single NWChem CDFT input deck for the pentacene example.

Reads the same YAML settings used by the trajectory runner and writes the
first selected frame to ``input/nwchem_cdft.nw``.
"""

from pathlib import Path

from q4mdkit.analysis.nwchem_cdft import (
    _extract_atom_types_from_topology,
    _iter_qm_coordinates,
    _load_qm_indices,
    load_config,
    render_nwchem_input,
)


CONFIG = Path(__file__).parent / "nwchem_cdft_pentacene.yaml"
OUTPUT = Path(__file__).parent / "input" / "nwchem_cdft.nw"


def main() -> None:
    config = load_config(CONFIG)
    qm_indices = _load_qm_indices(config.qm_atoms_file)
    atom_types = _extract_atom_types_from_topology(config.topology_path, qm_indices)

    frames = _iter_qm_coordinates(
        config.traj_path, config.topology_path, qm_indices
    )
    for frame_id, _time, coords in frames:
        if frame_id < config.start_frame:
            continue

        text = render_nwchem_input(
            atom_types=atom_types,
            coords_ang=coords,
            constraints=config.constraints,
            basis=config.basis,
            xc=config.xc,
            charge=config.charge,
            multiplicity=config.multiplicity,
            point_charges=None,
        )
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(text, encoding="utf-8")
        print(f"wrote {OUTPUT}")
        return

    raise SystemExit("no frames available")


if __name__ == "__main__":
    main()
