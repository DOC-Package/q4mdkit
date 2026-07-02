# TD-DFTB Simple Excitation Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove TD-DFTB orbital/state tracking and write up to `n_excitations` energy/oscillator-strength pairs horizontally for each successful frame.

**Architecture:** Keep the existing DFTB+ execution, EXC.DAT parsing, and optional spectrum calculation. Simplify the result path so raw excitation order flows directly from the parser to a single horizontal writer, with no eigenvector, overlap, transition, or tracking pipeline.

**Tech Stack:** Python 3, NumPy, PyYAML, pytest, DFTB+ HSD configuration

---

## File Structure

- `q4mdkit/analysis/tddftb.py`: retain TD-DFTB execution/parsing/spectrum behavior; remove tracking types and runtime branches; provide horizontal excitation output.
- `tests/test_tddftb_tracking.py`: replace tracking-only tests with output and configuration regression tests. Keeping the existing path avoids deleting a user-modified test file before its changes are incorporated.
- `examples/analysis/tddftb_settings.yaml`: remove tracked orbital settings from the tracked example configuration.
- `examples/quinacridone/analysis_tddftb/tddftb_settings.yaml`: verify it requires no removed keys; modify only if a removed key is present.

### Task 1: Specify Horizontal Excitation Output

**Files:**
- Modify: `tests/test_tddftb_tracking.py`
- Modify: `q4mdkit/analysis/tddftb.py`

- [ ] **Step 1: Replace tracking tests with a failing horizontal-output test**

```python
from q4mdkit.analysis.tddftb import (
    Excitation,
    append_excitation_frame,
    initialize_excitation_output,
)


def test_excitation_output_writes_energy_strength_pairs_horizontally(tmp_path):
    output_path = tmp_path / "excitations.dat"
    frame_data = {
        "frame": 5,
        "time_fs": 20.0,
        "excitations": [
            Excitation(energy_ev=1.1, oscillator_strength=0.2),
            Excitation(energy_ev=1.3, oscillator_strength=0.4),
        ],
    }

    initialize_excitation_output(output_path, n_excitations=2)
    append_excitation_frame(output_path, frame_data, n_excitations=2)

    assert output_path.read_text().splitlines() == [
        "# TD-DFTB Excitation Data",
        "# Frame  Time(fs)  E1(eV)  f1  E2(eV)  f2",
        "     5  20.00  1.1000  0.200000  1.3000  0.400000",
    ]
```

- [ ] **Step 2: Run the test and verify RED**

Run: `pytest -q tests/test_tddftb_tracking.py::test_excitation_output_writes_energy_strength_pairs_horizontally`

Expected: FAIL because the current writer takes `track_orbital_order` and emits one row per state with wavelength.

- [ ] **Step 3: Implement the minimal horizontal writer**

```python
def initialize_excitation_output(output_path: Path, n_excitations: int) -> None:
    with open(output_path, "w") as f:
        f.write("# TD-DFTB Excitation Data\n")
        columns = ["# Frame", "Time(fs)"]
        for state_index in range(1, n_excitations + 1):
            columns.extend([f"E{state_index}(eV)", f"f{state_index}"])
        f.write("  ".join(columns) + "\n")


def append_excitation_frame(
    output_path: Path,
    frame_data: Dict,
    n_excitations: int,
) -> None:
    fields = [f"{frame_data['frame']:6d}", f"{frame_data['time_fs']:.2f}"]
    for excitation in frame_data["excitations"][:n_excitations]:
        fields.extend([
            f"{excitation.energy_ev:.4f}",
            f"{excitation.oscillator_strength:.6f}",
        ])
    with open(output_path, "a") as f:
        f.write("  ".join(fields) + "\n")
```

- [ ] **Step 4: Run the test and verify GREEN**

Run: `pytest -q tests/test_tddftb_tracking.py::test_excitation_output_writes_energy_strength_pairs_horizontally`

Expected: `1 passed`.

- [ ] **Step 5: Add and verify truncation behavior**

```python
def test_excitation_output_limits_pairs_to_n_excitations(tmp_path):
    output_path = tmp_path / "excitations.dat"
    frame_data = {
        "frame": 0,
        "time_fs": 0.0,
        "excitations": [
            Excitation(energy_ev=1.0, oscillator_strength=0.1),
            Excitation(energy_ev=2.0, oscillator_strength=0.2),
            Excitation(energy_ev=3.0, oscillator_strength=0.3),
        ],
    }

    initialize_excitation_output(output_path, n_excitations=2)
    append_excitation_frame(output_path, frame_data, n_excitations=2)

    data_line = output_path.read_text().splitlines()[-1]
    assert data_line.split() == ["0", "0.00", "1.0000", "0.100000", "2.0000", "0.200000"]
```

Run: `pytest -q tests/test_tddftb_tracking.py`

Expected: `2 passed`.

### Task 2: Remove Tracking Configuration and Runtime

**Files:**
- Modify: `tests/test_tddftb_tracking.py`
- Modify: `q4mdkit/analysis/tddftb.py`

- [ ] **Step 1: Add a failing configuration regression test**

```python
from dataclasses import fields
from q4mdkit.analysis.tddftb import TDDFTBConfig


def test_tddftb_config_has_no_tracking_fields():
    field_names = {field.name for field in fields(TDDFTBConfig)}
    assert not any("tracking" in name or name.startswith("track_") for name in field_names)
    assert "state_timeseries_file" not in field_names
```

- [ ] **Step 2: Run the test and verify RED**

Run: `pytest -q tests/test_tddftb_tracking.py::test_tddftb_config_has_no_tracking_fields`

Expected: FAIL listing current tracking fields.

- [ ] **Step 3: Remove tracking-only production code**

In `q4mdkit/analysis/tddftb.py`:

- Remove imports of `read_eigenvectors_from_dir_spin_polarized` and `compute_cross_overlap_odin`.
- Remove `state_timeseries_file` and all tracking fields from `TDDFTBConfig` and `load_tddftb_config`.
- Reduce `Excitation` to energy, oscillator strength, optional wavelength, and raw dominant-transition fields used by EXC.DAT parsing.
- Remove `TransitionContribution`, TRA.DAT parsing/attachment, band/eigenvector readers, assignment helpers, tracking result types, trackers, overlap helpers, tracking application helpers, and tracked time-series writers.
- In `run_tddftb_analysis`, remove tracker creation and all frame-to-frame tracking state. Parse EXC.DAT directly and call:

```python
initialize_excitation_output(config.excitations_file, config.n_excitations)
append_excitation_frame(
    config.excitations_file,
    frame_exc_data,
    config.n_excitations,
)
```

- Keep wavelength calculation inside `Excitation` because optional spectrum output and the S1 summary still use the energy-to-wavelength conversion; do not write wavelength to `excitations.dat`.

- [ ] **Step 4: Run focused tests and static checks**

Run: `pytest -q tests/test_tddftb_tracking.py`

Expected: `3 passed`.

Run: `python -m py_compile q4mdkit/analysis/tddftb.py`

Expected: exit status 0 with no output.

Run: `rg -n "FrontierOrbital|NTOState|track_orbital|track_excited|state_timeseries|odin_overlap|TRA.DAT" q4mdkit/analysis/tddftb.py`

Expected: no matches.

### Task 3: Remove Obsolete Example Settings and Verify Behavior

**Files:**
- Modify: `examples/analysis/tddftb_settings.yaml`
- Modify if needed: `examples/quinacridone/analysis_tddftb/tddftb_settings.yaml`
- Do not modify: untracked `examples/pdi/analysis/` files

- [ ] **Step 1: Remove checked-in tracking settings**

Delete `track_orbital_order` and every `orbital_tracking_*`, `track_excited_state_order`, `excited_state_tracking_*`, and nested tracking-only `odin` block from tracked example YAML files. Keep `n_excitations`, broadening, energy range, spectrum settings, and unrelated output settings unchanged.

- [ ] **Step 2: Confirm no tracked example advertises removed settings**

Run: `git ls-files 'examples/**/*.yaml' | xargs rg -n "track_orbital|track_excited|orbital_tracking|excited_state_tracking|state_timeseries_file"`

Expected: no matches.

- [ ] **Step 3: Run the relevant and full test suites**

Run: `pytest -q tests/test_tddftb_tracking.py tests/test_cdftbci_overlap_dump.py tests/test_phase_tracking.py`

Expected: all selected tests pass.

Run: `pytest -q`

Expected: all tests pass. If unrelated environment-dependent tests fail, record exact failures and verify the TD-DFTB tests remain green.

- [ ] **Step 4: Review only intended diffs**

Run: `git diff --check`

Expected: no whitespace errors.

Run: `git diff -- q4mdkit/analysis/tddftb.py tests/test_tddftb_tracking.py examples/analysis/tddftb_settings.yaml examples/quinacridone/analysis_tddftb/tddftb_settings.yaml`

Expected: tracking removal, horizontal output, focused tests, and obsolete setting removal only; preserve pre-existing average-spectrum changes.
