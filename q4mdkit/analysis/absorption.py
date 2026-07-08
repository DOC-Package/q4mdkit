"""
Absorption spectrum from a spectral density J(omega) within the
second-order cumulant (Franck-Condon) approximation.

I(w)  proportional to  w * Re int_0^inf dt exp(i w t) chi(t)

chi(t) = sum_m |mu_m|^2 exp(-i eps_m t - g_m^(2)(t))

g_m^(2)(t) = (1/pi) int_0^inf dw J_m(w)/w^2
                 * [ coth(beta w / 2) (1 - cos w t) + i (sin w t - w t) ]

All quantities are expected in a single consistent unit system
(e.g. atomic units, where hbar = k_B = 1, energy and frequency
share the same unit). The user is responsible for supplying
``beta = 1/(k_B T)`` in inverse-energy units consistent with
``omega_J``, ``eps_m``, and the desired output frequency grid.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

__all__ = [
    "lineshape_function",
    "dipole_correlation",
    "absorption_spectrum",
    "oscillator_strength_to_mu2",
]

EV_PER_HARTREE = 27.211386245988


def oscillator_strength_to_mu2(oscillator_strength: float, energy_ev: float) -> float:
    """Convert dimensionless oscillator strength to |mu|^2 in atomic units."""
    if oscillator_strength < 0.0:
        raise ValueError("oscillator strength must be nonnegative")
    if energy_ev <= 0.0:
        raise ValueError("excitation energy must be positive")
    energy_hartree = energy_ev / EV_PER_HARTREE
    return 3.0 * oscillator_strength / (2.0 * energy_hartree)


def _as_2d_J(J: np.ndarray, n_modes: int, n_w: int) -> np.ndarray:
    """Broadcast J to shape (n_modes, n_w)."""
    J = np.asarray(J, dtype=float)
    if J.ndim == 1:
        if J.shape[0] != n_w:
            raise ValueError(
                f"J has length {J.shape[0]}, expected {n_w} (len(omega_J))."
            )
        return np.broadcast_to(J, (n_modes, n_w))
    if J.ndim == 2:
        if J.shape != (n_modes, n_w):
            raise ValueError(
                f"J has shape {J.shape}, expected ({n_modes}, {n_w})."
            )
        return J
    raise ValueError("J must be 1D (shared) or 2D (per-mode).")


def lineshape_function(
    t: np.ndarray,
    omega_J: np.ndarray,
    J: np.ndarray,
    beta: float,
    n_modes: int | None = None,
) -> np.ndarray:
    """
    Compute the second-order cumulant lineshape function g_m^(2)(t).

    Parameters
    ----------
    t : (Nt,) array
        Time grid (>= 0 typically).
    omega_J : (Nw,) array
        Evenly-spaced positive-frequency grid for the spectral density.
        omega_J[0] may be 0 (handled by removing the singular point).
    J : (Nw,) or (n_modes, Nw) array
        Spectral density J_m(omega). If 1D, the same J is used for every
        mode.
    beta : float
        Inverse temperature in units inverse to ``omega_J``.
    n_modes : int, optional
        Number of modes when ``J`` is 1D. Defaults to 1.

    Returns
    -------
    g : (n_modes, Nt) complex array
        Lineshape function for each mode.
    """
    t = np.asarray(t, dtype=float)
    w = np.asarray(omega_J, dtype=float)
    if w.ndim != 1 or w.size < 2:
        raise ValueError("omega_J must be a 1D array with at least 2 points.")

    dw = w[1] - w[0]
    if not np.allclose(np.diff(w), dw):
        raise ValueError("omega_J must be evenly spaced.")

    if n_modes is None:
        n_modes = 1 if np.asarray(J).ndim == 1 else np.asarray(J).shape[0]

    Jm = _as_2d_J(J, n_modes, w.size)

    # Drop omega = 0 to avoid the 1/w^2 singularity (J(0) = 0 physically).
    mask = w > 0.0
    w_pos = w[mask]
    Jm_pos = Jm[:, mask]

    coth = 1.0 / np.tanh(0.5 * beta * w_pos)  # (Nw,)
    pref = Jm_pos / (w_pos ** 2)

    # Bound temporary kernel memory for long trajectories by integrating
    # fixed-size blocks of time points.
    g = np.empty((n_modes, t.size), dtype=complex)
    time_chunk_size = 64
    for start in range(0, t.size, time_chunk_size):
        stop = min(start + time_chunk_size, t.size)
        wt = np.outer(t[start:stop], w_pos)
        real_kernel = coth * (1.0 - np.cos(wt))
        imag_kernel = np.sin(wt) - wt
        kernel = real_kernel + 1j * imag_kernel
        integrand = pref[:, None, :] * kernel[None, :, :]
        g[:, start:stop] = (
            np.trapezoid(integrand, w_pos, axis=-1) / np.pi
        )
    return g


def dipole_correlation(
    t: np.ndarray,
    mu2: Sequence[float],
    eps: Sequence[float],
    omega_J: np.ndarray,
    J: np.ndarray,
    beta: float,
) -> np.ndarray:
    """
    Compute chi(t) = sum_m |mu_m|^2 exp(-i eps_m t - g_m^(2)(t)).

    Parameters
    ----------
    t : (Nt,) array
        Time grid.
    mu2 : (n_modes,) array
        Squared transition dipole magnitudes |mu_m|^2.
    eps : (n_modes,) array
        Vertical excitation energies eps_m.
    omega_J, J, beta :
        See :func:`lineshape_function`.

    Returns
    -------
    chi : (Nt,) complex array
    """
    mu2 = np.asarray(mu2, dtype=float)
    eps = np.asarray(eps, dtype=float)
    if mu2.shape != eps.shape or mu2.ndim != 1:
        raise ValueError("mu2 and eps must be 1D arrays of equal length.")
    n_modes = mu2.size

    t = np.asarray(t, dtype=float)
    g = lineshape_function(t, omega_J, J, beta, n_modes=n_modes)  # (M, Nt)

    phase = np.exp(-1j * np.outer(eps, t) - g)        # (M, Nt)
    chi = np.einsum("m,mt->t", mu2, phase)            # (Nt,)
    return chi


def absorption_spectrum(
    mu2: Sequence[float],
    eps: Sequence[float],
    omega_J: np.ndarray,
    J: np.ndarray,
    beta: float,
    t_max: float,
    dt: float,
    damping: float = 0.0,
    omega_out: np.ndarray | None = None,
):
    """
    Numerically evaluate the absorption lineshape

        I(w)  =  w * Re int_0^{t_max} dt exp(i w t - damping*t) chi(t).

    Parameters
    ----------
    mu2, eps : (n_modes,) arrays
        Oscillator strengths |mu_m|^2 and vertical energies eps_m.
    omega_J : (Nw,) array
        Evenly-spaced frequency grid for J (same units as eps).
    J : (Nw,) or (n_modes, Nw) array
        Spectral density. 1D arrays are shared across modes.
    beta : float
        Inverse temperature, inverse-energy units consistent with eps.
    t_max : float
        Upper limit of the time integral. Choose t_max ~ several times
        the longest decay time of |chi(t)|.
    dt : float
        Time step used to sample chi(t). Must satisfy
        dt < pi / max(|eps|, max(omega_J)) (Nyquist for the relevant
        spectral content).
    damping : float, optional
        Phenomenological exponential damping (Lorentzian broadening
        with HWHM = damping). Useful when t_max is finite. Default 0.
    omega_out : array, optional
        Output frequency grid. Defaults to the FFT-native grid
        ``2*pi*np.fft.rfftfreq(Nt, dt)``.

    Returns
    -------
    omega : (Nout,) array
        Output angular-frequency grid.
    I : (Nout,) array
        Absorption lineshape I(w) (real, in arbitrary units).
    chi : (Nt,) complex array
        The dipole-dipole correlation function actually used.
    """
    if dt <= 0 or t_max <= 0:
        raise ValueError("dt and t_max must be positive.")

    Nt = int(np.ceil(t_max / dt)) + 1
    t = np.arange(Nt) * dt

    chi = dipole_correlation(t, mu2, eps, omega_J, J, beta)
    chi_d = chi * np.exp(-damping * t)

    # Trapezoidal endpoint correction for the half-line integral.
    chi_d = chi_d.copy()
    chi_d[0] *= 0.5
    chi_d[-1] *= 0.5

    if omega_out is None:
        # Default output grid: spans the relevant spectral range with
        # resolution ~ 2*pi / t_max.
        eps_arr = np.asarray(eps, dtype=float)
        w_span = float(np.max(np.abs(eps_arr))) + 10.0 * float(omega_J[-1])
        n_out = max(int(2.0 * w_span * Nt * dt / np.pi), 1024)
        omega = np.linspace(0.0, w_span, n_out)
    else:
        omega = np.asarray(omega_out, dtype=float)

    # Direct DFT on arbitrary output energies, chunked to bound memory use:
    #   int_0^T dt e^{i w t} f(t)  ~  dt * sum_n e^{i w t_n} f_n
    # phase has shape (Nout, Nt). For very large grids this is O(Nout*Nt)
    # but is conceptually simple and works on arbitrary omega grids.
    spec = np.empty(omega.size, dtype=complex)
    chunk_size = 256
    for start in range(0, omega.size, chunk_size):
        stop = min(start + chunk_size, omega.size)
        phase = np.exp(1j * np.outer(omega[start:stop], t))
        spec[start:stop] = phase @ chi_d * dt
    I = omega * np.real(spec)

    return omega, I, chi
