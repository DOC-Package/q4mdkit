#!/usr/bin/env python3
"""Plot J_lowdin (meV) vs Time (fs) from cdftbci_extracted.dat."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

here = Path(__file__).parent
data = np.loadtxt(here / "cdftbci_extracted.dat")
time = data[:, 1]
j_lowdin = data[:, 2]
mask = np.isfinite(j_lowdin)
time = time[mask]
j_lowdin = j_lowdin[mask]

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(time, j_lowdin, lw=0.8, color="tab:blue")
ax.set_xlabel("Time (fs)")
ax.set_ylabel(r"$J_\mathrm{lowdin}$ (meV)")
ax.set_title("CDFTB-CI electronic coupling vs time")
ax.grid(True, alpha=0.3)
ax.axhline(np.mean(j_lowdin), color="tab:red", ls="--", lw=0.8,
           label=f"mean = {np.mean(j_lowdin):.2f} meV")
ax.legend()
fig.tight_layout()

out = here / "jlowdin_vs_time.png"
fig.savefig(out, dpi=200)
print(f"Saved: {out}")
print(f"N={len(time)}, mean={np.mean(j_lowdin):.3f} meV, "
      f"std={np.std(j_lowdin):.3f} meV, "
      f"min={np.min(j_lowdin):.3f}, max={np.max(j_lowdin):.3f}")

# |J_lowdin| plot
abs_j = np.abs(j_lowdin)
fig2, ax2 = plt.subplots(figsize=(8, 4))
ax2.plot(time, abs_j, lw=0.8, color="tab:green")
ax2.set_xlabel("Time (fs)")
ax2.set_ylabel(r"$|J_\mathrm{lowdin}|$ (meV)")
ax2.set_title("CDFTB-CI |electronic coupling| vs time")
ax2.grid(True, alpha=0.3)
ax2.axhline(np.mean(abs_j), color="tab:red", ls="--", lw=0.8,
            label=f"mean = {np.mean(abs_j):.2f} meV")
ax2.legend()
fig2.tight_layout()

out2 = here / "jlowdin_abs_vs_time.png"
fig2.savefig(out2, dpi=200)
print(f"Saved: {out2}")
print(f"|J|: mean={np.mean(abs_j):.3f} meV, "
      f"std={np.std(abs_j):.3f} meV, "
      f"min={np.min(abs_j):.3f}, max={np.max(abs_j):.3f}")
