# TD-DFTB Simple Excitation Output Design

## Goal

Remove orbital and excited-state tracking from `q4mdkit.analysis.tddftb` for now. For each successful trajectory frame, write the configured number of raw TD-DFTB excitation energies and oscillator strengths in a compact horizontal format.

## Configuration

- `tddftb.n_excitations` remains the single setting that controls both DFTB+'s `NrOfExcitations` and the maximum number of excitation pairs written per frame.
- Remove all orbital-tracking and excited-state-tracking settings.
- Remove `output.state_timeseries_file` because there is no tracked time-series output.
- Preserve average-spectrum and individual-spectrum settings and behavior.

## Output

`output.excitations_file` contains one row per successful frame:

```text
# Frame  Time(fs)  E1(eV)  f1  E2(eV)  f2  ...  EN(eV)  fN
     0  0.00       1.2345  0.012345  1.5678  0.023456
```

Excitations retain the raw order returned by DFTB+, normally ascending excitation energy. Each energy is immediately followed by its oscillator strength. Wavelengths, transition labels, tracking diagnostics, and tracked-state indices are not written.

The writer outputs at most `n_excitations` pairs. If DFTB+ returns fewer states, it writes the available pairs without padding.

## Implementation Boundaries

- Keep `Excitation`, EXC.DAT parsing, spectrum generation, HSD preparation, trajectory iteration, and summary statistics.
- Remove tracking-only imports, data fields, parsers, classes, helper functions, runtime state, and output functions.
- Replace the current vertical excitation writer with a horizontal writer that accepts `n_excitations`.
- Remove or replace tracking-specific tests with focused tests for horizontal output and truncation.
- Update checked-in example settings that expose removed tracking keys.

## Error Handling

Existing behavior remains: failed DFTB+ frames and frames without parsed excitations are skipped. The excitation output is initialized before frame processing and remains valid even when no frame succeeds.

## Verification

- A unit test proves one frame is written as `E1, f1, E2, f2` on one line.
- A unit test proves output is truncated to `n_excitations`.
- Configuration tests or direct assertions confirm removed tracking settings are no longer part of `TDDFTBConfig`.
- Run the targeted TD-DFTB tests, then the full test suite.
