#!/usr/bin/env python3
"""Generate a REAL energy-conservation figure from transient_energy_2d
(replaces the fabricated np.ones Fig. 3). Shows E_H(t)/E_H(0) and its
machine-precision drift for the symmetric Ventcel leapfrog scheme.

Output: ../figures/fig3_energy.pdf and .png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _d in ('solvers', 'utils'):
    _p = _os.path.join(_SRC, _d)
    if _p not in _sys.path: _sys.path.insert(0, _p)
del _sys, _os, _SRC, _d, _p
import pubstyle; pubstyle.apply()
import matplotlib.pyplot as plt
from transient_energy_2d import run


def make_figure():
    drift, E, E0 = run(refine=5, n_steps=4000, report=True)
    rel = E / E0
    t = np.arange(len(E))

    fig, ax = plt.subplots(2, 1, figsize=(3.8, 5.4))

    # panel 0: discrete energy ratio  (single curve -> BW[0])
    # dense 4000-step trace: no markers (explicit "none" overrides the BW cycler)
    s0 = dict(pubstyle.BW[0]); s0["marker"] = "none"
    ax[0].plot(t, rel, lw=0.9, **s0)
    ax[0].set_xlabel("time step $n$")
    ax[0].set_ylabel(r"$E_H(t)/E_H(0)$")
    ax[0].set_title("Discrete energy (symmetric Ventcel leapfrog)")
    ax[0].set_ylim(1 - 2e-14, 1 + 2e-14)
    ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # panel 1: drift on log scale  (single curve -> BW[1])
    s1 = dict(pubstyle.BW[1]); s1["marker"] = "none"
    ax[1].semilogy(t, np.abs(rel - 1.0) + 1e-18, lw=0.9, **s1)
    ax[1].set_xlabel("time step $n$")
    ax[1].set_ylabel(r"$|E_H(t)/E_H(0)-1|$")
    ax[1].set_title(f"drift (max = {drift:.1e}, machine precision)")
    ax[1].axhline(1e-8, color="0.0", ls=":", lw=0.8)
    fig.tight_layout()

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
    os.makedirs(out, exist_ok=True)
    base = os.path.join(out, "fig3_energy")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved REAL energy figure: {base}.pdf / .png")
    print(f"  (max drift over 4000 steps = {drift:.3e}; genuine roundoff, not np.ones)")


if __name__ == "__main__":
    make_figure()
