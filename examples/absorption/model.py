"""
Model absorption spectrum within the second-order cumulant approximation.

Demonstrates :func:`q4mdkit.analysis.absorption.absorption_spectrum`
on a simple toy system:

  * single electronic state (M = 1)
  * Ohmic spectral density with exponential cutoff
        J(w) = (pi/2) * lambda * (w / w_c) * exp(-w / w_c)
    so that the reorganization energy is ``lambda``.
  * physical units: electron volts (eV) for energy/frequency,
    femtoseconds (fs) for time. We work in a unit system where
    hbar = 1 and 1 eV * 1 fs = hbar / (hbar / (eV*fs)) = 0.6582 (=> 1/hbar_eVfs).
    To keep the formulas dimensionally simple we use:
        time:        fs
        energy:      eV
        with omega in (1/fs) = energy_in_eV / hbar_eVfs
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from q4mdkit.analysis.absorption import absorption_spectrum


# ---------------------------------------------------------------------------
# Unit conventions
# ---------------------------------------------------------------------------
HBAR_EV_FS = 0.6582119569       # hbar in eV * fs
KB_EV_PER_K = 8.617333262e-5    # Boltzmann constant, eV / K


def ev_to_invfs(e_ev: np.ndarray | float) -> np.ndarray | float:
    """Convert energy in eV to angular frequency in 1/fs (omega = E/hbar)."""
    return np.asarray(e_ev) / HBAR_EV_FS


def invfs_to_ev(w_invfs: np.ndarray | float) -> np.ndarray | float:
    return np.asarray(w_invfs) * HBAR_EV_FS


# ---------------------------------------------------------------------------
# Model definition
# ---------------------------------------------------------------------------
def ohmic_spectral_density(
    omega: np.ndarray, reorg_energy: float, omega_c: float
) -> np.ndarray:
    """
    Ohmic J(w) with exponential cutoff. Reorganization energy

        lambda = (1/pi) * int_0^inf dw J(w) / w  =  reorg_energy.
    """
    J = np.zeros_like(omega)
    pos = omega > 0
    w = omega[pos]
    J[pos] = 0.5 * np.pi * reorg_energy * (w / omega_c) * np.exp(-w / omega_c)
    return J


def main() -> None:
    # -- Electronic structure (single mode) ---------------------------------
    eps_eV = np.array([2.50])              # vertical excitation energy
    mu2 = np.array([1.0])                  # |mu|^2 (arb. units)

    # -- Bath / spectral density --------------------------------------------
    reorg_eV = 0.10                        # reorganization energy
    omega_c_eV = 0.1                      # cutoff (~ 400 cm^-1)
    T_K = 10.0

    # convert to angular-frequency units (1/fs) for the time integral
    eps = ev_to_invfs(eps_eV)
    reorg = ev_to_invfs(reorg_eV)
    omega_c = ev_to_invfs(omega_c_eV)
    beta = HBAR_EV_FS / (KB_EV_PER_K * T_K)   # 1 / (k_B T) in fs

    # spectral-density grid (evenly spaced)
    n_w = 4001
    w_max = 20.0 * omega_c
    omega_J = np.linspace(0.0, w_max, n_w)
    J = ohmic_spectral_density(omega_J, reorg, omega_c)

    # -- Time-domain integral parameters ------------------------------------
    dt = 0.10                              # fs
    t_max = 2000.0                         # fs
    damping = 1.0 / 50.0                  # 1/fs  (-> ~1.3 meV HWHM)

    # -- Output frequency grid ----------------------------------------------
    omega_out_eV = np.linspace(1.5, 3.5, 2001)
    omega_out = ev_to_invfs(omega_out_eV)

    omega, I, chi = absorption_spectrum(
        mu2=mu2,
        eps=eps,
        omega_J=omega_J,
        J=J,
        beta=beta,
        t_max=t_max,
        dt=dt,
        damping=damping,
        omega_out=omega_out,
    )

    # Normalize for plotting
    I_plot = I / np.max(np.abs(I))

    # -- Plot ---------------------------------------------------------------
    t = np.arange(int(np.ceil(t_max / dt)) + 1) * dt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    ax.plot(t, np.real(chi), label="Re chi(t)")
    ax.plot(t, np.imag(chi), label="Im chi(t)", linestyle="--")
    ax.plot(t, np.abs(chi), label="|chi(t)|", color="k", linewidth=0.8)
    ax.set_xlabel("t  [fs]")
    ax.set_ylabel("chi(t)")
    ax.set_xlim(0, min(t_max, 1500.0))
    ax.legend()
    ax.set_title("Dipole autocorrelation")

    ax = axes[1]
    ax.plot(omega_out_eV, I_plot)
    ax.axvline(eps_eV[0], color="gray", linestyle=":", label="vertical eps")
    ax.axvline(eps_eV[0] - reorg_eV, color="red", linestyle=":",
               label="0-0 (eps - lambda)")
    ax.set_xlabel("omega  [eV]")
    ax.set_ylabel("I(omega) (norm.)")
    ax.set_title(
        f"Absorption (T={T_K:.0f} K, lambda={reorg_eV:.2f} eV, "
        f"omega_c={omega_c_eV:.2f} eV)"
    )
    ax.legend()

    fig.tight_layout()
    out_png = "absorption_model.png"
    fig.savefig(out_png, dpi=150)
    print(f"Saved: {out_png}")

    # also dump the spectrum
    np.savetxt(
        "absorption_model.dat",
        np.column_stack([omega_out_eV, I]),
        header="omega[eV]  I(omega)[arb]",
    )
    print("Saved: absorption_model.dat")


if __name__ == "__main__":
    main()
