#!/usr/bin/env python3
"""
Animation of the resonant band gap (negative effective mass mechanism).

The animation sweeps frequency ω from 0 to 1.6 ω₁ and shows:
  Left  — effective surface mass M_eff(ω) being traced; gap region in red
  Right — wave snapshot U(x, t=fixed) from a 1D toy model: propagating
          below ω₁, evanescent (decaying) inside the gap, propagating again above

Parameters from resonant_mass_2d.py (soft dense inclusion, locally-resonant).
Output: figures/anim_bandgap.gif   (~120 frames, 6 s at 20 fps)
"""
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _d in ('solvers', 'utils'):
    _p = _os.path.join(_SRC, _d)
    if _p not in _sys.path: _sys.path.insert(0, _p)
del _sys, _os, _SRC, _d, _p

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import pubstyle; pubstyle.apply()
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
import resonant_mass_2d as rm

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)

# ── exact parameters from resonant_mass_2d (same as fig_bandgap.pdf, Fig. 9) ─
MU_M  = rm.MU_M
RHO_M = rm.RHO_M

print("  Computing inclusion modes (FE eigenproblem)...")
_omega_n, _a_n, _M_incl, _ = rm.inclusion_modes()
omega1 = _omega_n[0]   # fundamental resonance (FE, verified vs Bessel j_{0,1})

# Effective surface mass — same formula and same data as make_figs_modelA::fig_bandgap()
def M_eff(omega):
    with np.errstate(divide='ignore', invalid='ignore'):
        result = rm.effective_mass(np.atleast_1d(np.real(omega)), _omega_n, _a_n, _M_incl)
        return result + 0j     # always an array, even for scalar omega

# 1-D wave field model: u(x) = exp(i k x), k² = rho_m ω² / mu_m - M_eff ω² / (mu_m * h)
# Simplified: in the gap M_eff < 0 => k imaginary => evanescent
def wave_snapshot(omega, x, h=0.01):
    Me = float(np.real(M_eff(omega)).ravel()[0])
    k2 = float((RHO_M * omega**2 - Me * omega**2 / h) / MU_M)
    if k2 > 0:
        k = np.sqrt(k2)
        return np.real(np.exp(1j * k * x))   # propagating
    elif k2 < 0:
        k = np.sqrt(-k2)
        return np.exp(-k * x)                # evanescent
    else:
        return np.ones_like(x)

# ── precompute full M_eff curve ───────────────────────────────────────────────
OMEGA_MAX = 1.6 * omega1
# avoid singularity at omega1
omega_lo = np.linspace(1e-3 * omega1, 0.97 * omega1, 400)
omega_hi = np.linspace(1.03 * omega1, OMEGA_MAX, 200)
omega_all = np.concatenate([omega_lo, omega_hi])
Me_all = np.real(M_eff(omega_all))

# clip for display
MCLIP = 3 * abs(_M_incl)
Me_clipped = np.clip(Me_all, -MCLIP, MCLIP)

# gap extent
gap_lo = omega1
gap_hi_idx = np.where((omega_hi > omega1) & (np.real(M_eff(omega_hi)) > 0))[0]
gap_hi = omega_hi[gap_hi_idx[0]] if gap_hi_idx.size else OMEGA_MAX

N_FRAMES = 120
omega_anim = np.concatenate([
    np.linspace(1e-3 * omega1, 0.95 * omega1, 40),
    np.linspace(1.05 * omega1, OMEGA_MAX, 30),
])
x_wave = np.linspace(0, 1, 300)


def make_animation():
    print("  Building band-gap animation...")

    fig, (ax_m, ax_w) = plt.subplots(1, 2, figsize=(4.8, 2.2))

    # ── left: M_eff(ω) ──────────────────────────────────────────────────────
    freq_khz = omega_all / (2 * np.pi * 1e3)
    f1_khz   = omega1   / (2 * np.pi * 1e3)
    fmax_khz = OMEGA_MAX / (2 * np.pi * 1e3)

    # shaded gap region
    ax_m.axvspan(f1_khz, gap_hi / (2 * np.pi * 1e3),
                 alpha=0.15, color="red", zorder=0, label="band gap")
    ax_m.axhline(0, color="0.6", lw=0.6, ls=":")
    ax_m.axvline(f1_khz, color="0.4", lw=0.7, ls="--", label=r"$\omega_1$")

    # full curve (faint grey, drawn once)
    for seg_omega, seg_Me in [(omega_lo, Me_clipped[:len(omega_lo)]),
                               (omega_hi, Me_clipped[len(omega_lo):])]:
        ax_m.plot(seg_omega / (2 * np.pi * 1e3), seg_Me,
                  color="0.80", lw=0.8, zorder=1)

    line_m,  = ax_m.plot([], [], color=pubstyle.BW[0]["color"] if hasattr(pubstyle.BW[0], '__getitem__') else "k",
                          lw=1.2, zorder=3)
    dot_m,   = ax_m.plot([], [], "o", ms=5, color="tab:red", zorder=4)

    ax_m.set_xlim(0, fmax_khz)
    ax_m.set_ylim(-MCLIP * 1.1, MCLIP * 1.1)
    ax_m.set_xlabel(r"$f$ (kHz)")
    ax_m.set_ylabel(r"$M_\mathrm{eff}(\omega)$ (kg/m²)")
    ax_m.set_title("Effective surface mass")
    ax_m.legend(fontsize=6, loc="upper right")

    # annotate
    ax_m.text(f1_khz * 0.5, MCLIP * 0.7, "propagating\n$M_\\mathrm{eff}>0$",
              ha="center", va="center", fontsize=5.5, color="0.3")
    ax_m.text(f1_khz * 1.08, -MCLIP * 0.6, "gap\n$M_\\mathrm{eff}<0$",
              ha="left", va="center", fontsize=5.5, color="tab:red")

    # ── right: wave field ────────────────────────────────────────────────────
    line_w,  = ax_w.plot(x_wave, wave_snapshot(omega_anim[0], x_wave),
                         lw=1.0, **{k: v for k, v in pubstyle.BW[0].items() if k != "marker"})
    ax_w.axvline(0.0, color="k", lw=0.8, ls="--", label=r"$\Gamma$")
    ax_w.set_xlim(0, 1)
    ax_w.set_ylim(-1.2, 1.2)
    ax_w.set_xlabel(r"$x_1$")
    ax_w.set_ylabel(r"$U(x_1)$")
    ax_w.set_title("Wave field (1-D homogenized model)")
    ax_w.legend(fontsize=6)

    freq_txt  = ax_w.text(0.97, 0.97, "", transform=ax_w.transAxes,
                           ha="right", va="top", fontsize=6.5)
    state_txt = ax_w.text(0.97, 0.88, "", transform=ax_w.transAxes,
                           ha="right", va="top", fontsize=6.5)

    fig.tight_layout()

    # colour-map for M_eff line: blue=propagating, red=gap
    def _me_scalar(omega):
        """M_eff as a Python float for a single frequency."""
        return float(np.real(M_eff(omega)).ravel()[0])

    def frame_color(omega):
        return "tab:red" if _me_scalar(omega) < 0 else "tab:blue"

    def update(i):
        omega = omega_anim[i % len(omega_anim)]
        f_khz = omega / (2 * np.pi * 1e3)

        # trace M_eff up to this frame
        traced = omega_anim[:i + 1]
        Me_tr  = np.real(M_eff(traced))
        Me_tr  = np.clip(Me_tr, -MCLIP, MCLIP)
        line_m.set_data(traced / (2 * np.pi * 1e3), Me_tr)
        dot_m.set_data([f_khz], [np.clip(_me_scalar(omega), -MCLIP, MCLIP)])
        dot_m.set_color(frame_color(omega))

        # wave snapshot
        U = wave_snapshot(omega, x_wave)
        line_w.set_ydata(U)

        Me_val = np.real(M_eff(omega))
        in_gap = Me_val < 0
        freq_txt.set_text(f"$f$ = {f_khz:.1f} kHz")
        state_txt.set_text("EVANESCENT (gap)" if in_gap else "propagating")
        state_txt.set_color("tab:red" if in_gap else "tab:blue")

        # wave background colour
        ax_w.set_facecolor("#fff0f0" if in_gap else "#f0fff0")

        return line_m, dot_m, line_w, freq_txt, state_txt

    ani = animation.FuncAnimation(fig, update, frames=len(omega_anim),
                                  interval=50, blit=False)
    out = os.path.join(FIG, "anim_bandgap.gif")
    ani.save(out, writer=animation.PillowWriter(fps=20), dpi=60)
    plt.close(fig)
    print(f"  saved: {out}")


if __name__ == "__main__":
    make_animation()
