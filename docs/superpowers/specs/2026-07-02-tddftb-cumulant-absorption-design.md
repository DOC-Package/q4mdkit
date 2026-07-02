# TD-DFTB Energy-Gap Spectral Density and Cumulant Absorption Design

## Goal

Add two separate, reusable workflows:

1. Build a spectral density from a user-selected excited-state energy trajectory in horizontal `excitations.dat` output.
2. Calculate a second-order-cumulant absorption spectrum from that spectral density using the existing theory implementation in `q4mdkit.analysis.absorption`.

The default selected state is S5, but every positive one-based state index is supported.

## Components

### Spectral-density module

Create `q4mdkit/analysis/energy_gap_spectral_density.py` with pure functions to:

- parse `Frame Time(fs) E1 f1 ... EN fN` files by header labels;
- select `E<state>` and `f<state>` for a requested state;
- validate finite, uniformly spaced time samples and available state columns;
- subtract the mean gap to obtain `delta_E(t)`;
- construct the classical two-sided energy-gap power spectrum;
- apply the harmonic-bath classical quantum correction

  `J(omega) = beta * omega * S_deltaE(omega) / 2`;

- write the spectral density and metadata needed by the absorption calculation.

The PSD uses a mean-square-normalized Hann window and a real FFT. Internal natural units are energy in eV and time in `eV^-1`, using `hbar = 0.6582119569 eV fs`. For a trajectory sampled in fs, `dt_natural = dt_fs / hbar`. The FFT estimator is

`S_deltaE(omega_k) = dt_natural / (N * mean(window^2)) * abs(FFT(window * delta_E))^2`.

Only nonnegative FFT frequencies are written. The zero-frequency value of `J` is explicitly zero.

### Spectral-density CLI

Create `q4mdkit/analysis/build_excitation_spectral_density.py`.

Required/important options:

- positional input `excitations.dat`;
- `--state`, default `5`;
- `--temperature`, default `300 K`;
- `--output`, default `spectral_density_S<state>.dat`;
- optional zero-padding length for smoother frequency interpolation without claiming improved physical resolution.

Output columns:

- `omega_eV`;
- `wavenumber_cm-1`;
- `J_eV`.

Header metadata includes state index, temperature, sample count, `dt_fs`, mean excitation energy, mean oscillator strength, and energy standard deviation.

### Absorption CLI

Create `q4mdkit/analysis/calculate_cumulant_absorption.py`. It reads the spectral-density file and metadata, then calls the existing functions in `q4mdkit.analysis.absorption`.

It converts mean oscillator strength to the squared transition dipole in atomic units:

`mu2_au = 3 * mean_f / (2 * mean_energy_hartree)`.

The absolute spectrum remains in arbitrary units. The CLI exposes:

- positional spectral-density file;
- `--temperature`, defaulting to the file metadata;
- `--t-max-fs` and `--dt-fs` for the correlation-time integral;
- `--damping-ev`, default `0`, for optional phenomenological damping;
- `--energy-min`, `--energy-max`, and `--n-energy` for the output grid;
- `--output` for the absorption spectrum;
- optional `--correlation-output` for `chi(t)`.

Output columns:

- `energy_eV`;
- `wavelength_nm` (`inf` at zero energy);
- raw `intensity`;
- `normalized_intensity`, scaled by the maximum positive raw intensity.

### Existing absorption module

Retain the separation between construction and absorption. `q4mdkit/analysis/absorption.py` remains responsible only for `g^(2)(t)`, `chi(t)`, and `I(omega)`. Add validation or memory-safe chunking there only where required by tests and realistic CLI grids.

## Data Flow

```text
excitations.dat
  -> select E[state], f[state]
  -> delta_E(t)
  -> S_deltaE(omega)
  -> J_state(omega)
  -> spectral_density_Sn.dat
  -> g_state^(2)(t)
  -> chi(t)
  -> I(omega)
  -> absorption_Sn.dat
```

## Validation and Error Handling

- Reject state indices below one.
- Reject missing or incomplete `E<state>`/`f<state>` columns.
- Reject fewer than four samples, nonfinite values, nonmonotonic time, or nonuniform sampling.
- Reject nonpositive temperature and invalid frequency/time/output grids.
- Require spectral-density metadata for mean energy and oscillator strength unless explicit CLI overrides are supplied.
- Do not silently select the brightest state when a fixed state was requested.

## Tests

- Parse multiple selectable state columns from a synthetic horizontal excitation file.
- Reject unavailable state indices and irregular time grids.
- Verify a sinusoidal energy-gap fluctuation produces a nonnegative spectral-density peak at the expected frequency.
- Verify `J(0) = 0` and output metadata round-trips.
- Verify zero spectral density reduces the cumulant result to the expected undamped phase correlation.
- Verify oscillator-strength-to-dipole conversion.
- Exercise both CLIs on temporary files.
- Run the existing TD-DFTB and absorption-related tests plus the repository suite, excluding the three known legacy tests that require absent `.npy` fixtures.
