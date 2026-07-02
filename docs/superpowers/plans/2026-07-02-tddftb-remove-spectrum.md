# TD-DFTB Spectrum Removal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove every broadened-spectrum feature from TD-DFTB analysis while preserving horizontal excitation-energy and oscillator-strength output.

**Architecture:** Simplify `run_tddftb_analysis` to a direct trajectory → DFTB+ → EXC.DAT → `excitations.dat` pipeline. Delete spectrum-only configuration, data structures, numerical processing, output, and documentation rather than retaining disabled branches.

**Tech Stack:** Python 3, PyYAML, pytest, DFTB+

---

## File Structure

- `q4mdkit/analysis/tddftb.py`: remove spectrum configuration, Gaussian processing, wavelength conversion, and spectrum runtime branches.
- `tests/test_tddftb_tracking.py`: extend the existing configuration regression test to reject spectrum-only fields.
- `examples/analysis/tddftb_settings.yaml`: remove spectrum-only YAML settings.
- `examples/quinacridone/analysis_tddftb/tddftb_settings.yaml`: remove spectrum-only YAML settings.
- `examples/analysis/run_tddftb.py` and `examples/quinacridone/analysis_tddftb/run_tddftb.py`: remove average-spectrum claims from usage documentation.

### Task 1: Remove Spectrum API and Runtime

**Files:**
- Modify: `tests/test_tddftb_tracking.py`
- Modify: `q4mdkit/analysis/tddftb.py`

- [ ] **Step 1: Add a failing spectrum-field regression assertion**

Extend the configuration test with:

```python
    spectrum_fields = {
        "spectra_file",
        "output_average_spectrum",
        "output_individual_spectra",
        "broadening_ev",
        "energy_range_ev",
        "n_energy_points",
    }
    assert field_names.isdisjoint(spectrum_fields)
```

- [ ] **Step 2: Run the test and verify RED**

Run: `pytest -q tests/test_tddftb_tracking.py::test_tddftb_config_has_no_tracking_fields`

Expected: FAIL because all six spectrum-only fields still exist.

- [ ] **Step 3: Remove the spectrum implementation**

In `q4mdkit/analysis/tddftb.py`:

- Delete `NM_PER_EV`, `spectra_file`, `output_average_spectrum`, `output_individual_spectra`, `broadening_ev`, `energy_range_ev`, and `n_energy_points`.
- Delete parsing and constructor arguments for the corresponding YAML keys.
- Reduce `Excitation` to raw energy, oscillator strength, and optional transition metadata:

```python
@dataclass
class Excitation:
    energy_ev: float
    oscillator_strength: float
    transition_from: Optional[int] = None
    transition_to: Optional[int] = None
    transition_weight: Optional[float] = None
    ks_energy_ev: Optional[float] = None
```

- Delete `compute_absorption_spectrum`.
- In `run_tddftb_analysis`, remove energy/wavelength grids, `all_spectra`, Gaussian processing, individual/average file writes, peak output, and wavelength text from the S1 summary.
- Update the module and function docstrings to describe excitation extraction rather than absorption-spectrum generation.

- [ ] **Step 4: Run focused tests and static checks**

Run: `pytest -q tests/test_tddftb_tracking.py`

Expected: `3 passed`.

Run: `python -m py_compile q4mdkit/analysis/tddftb.py`

Expected: exit status 0.

Run: `rg -n "spectrum|spectra|broadening|wavelength|NM_PER_EV|all_spectra|avg_spectrum" q4mdkit/analysis/tddftb.py`

Expected: no matches.

### Task 2: Remove Obsolete Settings and Verify

**Files:**
- Modify: `examples/analysis/tddftb_settings.yaml`
- Modify: `examples/quinacridone/analysis_tddftb/tddftb_settings.yaml`
- Modify: `examples/analysis/run_tddftb.py`
- Modify: `examples/quinacridone/analysis_tddftb/run_tddftb.py`

- [ ] **Step 1: Remove spectrum-only example configuration and documentation**

Delete `spectra_file`, `average_spectrum`, `individual_spectra`, `broadening_ev`, `energy_range_ev`, and `n_energy_points` from both tracked TD-DFTB YAML examples. Remove average/individual absorption-spectrum output descriptions from the two runner docstrings while retaining the `excitations.dat` description.

- [ ] **Step 2: Confirm obsolete symbols are absent from tracked TD-DFTB files**

Run:

```bash
rg -n "spectra_file|average_spectrum|individual_spectra|broadening_ev|energy_range_ev|n_energy_points|compute_absorption_spectrum" \
  q4mdkit/analysis/tddftb.py \
  tests/test_tddftb_tracking.py \
  examples/analysis/tddftb_settings.yaml \
  examples/analysis/run_tddftb.py \
  examples/quinacridone/analysis_tddftb/tddftb_settings.yaml \
  examples/quinacridone/analysis_tddftb/run_tddftb.py
```

Expected: only the test's `spectrum_fields` names match; production and example files have no matches.

- [ ] **Step 3: Run verification**

Run: `pytest -q tests/test_tddftb_tracking.py tests/test_cdftbci_overlap_dump.py tests/test_phase_tracking.py`

Expected: `18 passed`.

Run: `pytest -q --ignore=examples/pentacene775-single/analysis/old/test_spectral_density_conversion.py --ignore=examples/pentacene775-single/analysis/old/test_temperature_prefactor.py --ignore=examples/pentacene775-single/analysis/old/test_window_effect.py`

Expected: all collected tests pass.

Run: `git diff --check`

Expected: no whitespace errors.
