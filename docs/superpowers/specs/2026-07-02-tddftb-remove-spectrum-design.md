# TD-DFTB Spectrum Removal Design

## Goal

Remove all broadened-spectrum generation from `q4mdkit.analysis.tddftb`. The analysis will only write the configured number of raw excitation energies and oscillator strengths for each successful trajectory frame.

This specification supersedes the spectrum-preservation requirements in `2026-07-02-tddftb-simple-excitation-output-design.md`.

## Configuration

Remove all spectrum-only configuration:

- `output.spectra_file`
- `output.average_spectrum`
- `output.individual_spectra`
- `tddftb.broadening_ev`
- `tddftb.energy_range_ev`
- `tddftb.n_energy_points`

The remaining TD-DFTB setting `n_excitations` controls both DFTB+'s requested roots and the number of energy/oscillator-strength pairs written.

## Runtime Behavior

- Do not construct an energy or wavelength grid.
- Do not Gaussian-broaden individual transitions.
- Do not retain per-frame spectra in memory.
- Do not write individual or average spectrum files.
- Do not calculate spectrum standard deviations or print an average-spectrum peak.
- Continue writing `excitations.dat` incrementally as `Frame Time E1 Osc1 ... EN OscN`.

## Code Removal

- Remove `compute_absorption_spectrum`.
- Remove spectrum-only fields and parsing from `TDDFTBConfig` and `load_tddftb_config`.
- Remove wavelength conversion from `Excitation` and the final summary because wavelength is not part of the requested output.
- Update module documentation and checked-in TD-DFTB examples so they no longer advertise spectrum output.
- Preserve EXC.DAT parsing, DFTB+ execution, trajectory iteration, and excitation output.

## Testing

- Add a regression assertion that `TDDFTBConfig` exposes no spectrum-only fields.
- Retain horizontal excitation-output and truncation tests.
- Verify no spectrum-generation symbols or obsolete YAML keys remain in tracked TD-DFTB files.
- Run targeted tests and the repository suite with the three previously identified environment-dependent legacy tests excluded.
