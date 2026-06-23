#!/usr/bin/env python3
"""
Animation of the 2D transient wave propagation (homogenized MsFEM medium).

Two panels:
  Left  — field U(x, y, t)  : Gaussian wave packet crossing interface Gamma={x=0.5}
  Right — discrete energy E(t)/E(0) accumulated in real time

Parameters identical to the paper's leapfrog solver (MU*, RHO* from cell problem).
Output: figures/anim_wave.gif   (~100 frames, ~5 s at 20 fps)
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
import matplotlib.tri as mtri
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm
from skfem.helpers import dot, grad

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)

# ── physical parameters (homogenized, from cell_problem_2d) ──────────────────
MU  = 1.747e10      # Pa
RHO = 3825.0        # kg/m³
C   = np.sqrt(MU / RHO)

# ── simulation parameters ────────────────────────────────────────────────────
REFINE     = 3      # 81 interior DOFs — fast and GIF-friendly
N_STEPS    = 360
SAVE_EVERY = 6      # → 60 frames
CFL        = 0.9


@BilinearForm
def mass_form(u, v, w):
    return RHO * u * v


@BilinearForm
def stiff_form(u, v, w):
    return MU * dot(grad(u), grad(v))


def solve():
    m = MeshTri().refined(REFINE)
    vb = Basis(m, ElementTriP1())
    M = mass_form.assemble(vb)
    K = stiff_form.assemble(vb)
    interior = vb.complement_dofs(vb.get_dofs())
    Mi = M[interior][:, interior].tocsc()
    Ki = K[interior][:, interior].tocsc()
    ml = np.asarray(Mi.sum(axis=1)).ravel()
    inv_ml = 1.0 / ml
    lam_max = spla.eigsh(Ki, k=1, M=sp.diags(ml), which="LA",
                         return_eigenvectors=False)[0]
    dt = CFL * 2.0 / np.sqrt(lam_max)

    # Gaussian wave packet centred at (0.25, 0.5), moving right
    px, py = m.p[0], m.p[1]
    sig = 0.07
    U0_full = np.exp(-((px - 0.25)**2 + (py - 0.5)**2) / (2 * sig**2))
    U0 = U0_full[interior]

    # leapfrog initialisation
    U_prev = U0.copy()
    acc0   = inv_ml * (-(Ki @ U0))
    U_curr = U0 + 0.5 * dt**2 * acc0

    def energy(Un, Unp1):
        V = (Unp1 - Un) / dt
        return 0.5 * np.dot(V, ml * V) + 0.5 * np.dot(Unp1, Ki @ Un)

    E0 = energy(U_prev, U_curr)
    snaps, times, E_hist = [], [], []

    for n in range(N_STEPS):
        acc    = inv_ml * (-(Ki @ U_curr))
        U_next = 2.0 * U_curr - U_prev + dt**2 * acc
        E_n    = energy(U_curr, U_next)
        E_hist.append(E_n / E0)
        if n % SAVE_EVERY == 0:
            U_full = np.zeros(m.p.shape[1])
            U_full[interior] = U_curr
            snaps.append(U_full.copy())
            times.append(n * dt)
        U_prev, U_curr = U_curr, U_next

    return m, snaps, times, np.array(E_hist), dt, E0


def make_animation():
    print("  Running transient solve for wave animation...")
    m, snaps, times, E_hist, dt, E0 = solve()

    tri = mtri.Triangulation(m.p[0], m.p[1], m.t.T)
    vlim = max(abs(s).max() for s in snaps) * 1.05
    n_frames = len(snaps)
    t_axis = np.arange(len(E_hist)) * dt * SAVE_EVERY

    fig, (ax_field, ax_en) = plt.subplots(1, 2, figsize=(4.8, 2.2))

    # ── left: wave field ────────────────────────────────────────────────────
    tcf = ax_field.tripcolor(tri, snaps[0], cmap=pubstyle.FIELD_CMAP,
                             vmin=-vlim, vmax=vlim, shading="gouraud")
    ax_field.axvline(0.5, color="k", lw=0.8, ls="--", label=r"$\Gamma$")
    ax_field.set_xlim(0, 1); ax_field.set_ylim(0, 1)
    ax_field.set_aspect("equal")
    ax_field.set_xlabel(r"$x_1$"); ax_field.set_ylabel(r"$x_2$")
    ax_field.set_title(r"$U(x,t)$ — homogenized MsFEM")
    ax_field.legend(loc="upper right", fontsize=6)
    fig.colorbar(tcf, ax=ax_field, fraction=0.046, pad=0.04)

    time_txt = ax_field.text(0.02, 0.97, "", transform=ax_field.transAxes,
                              va="top", fontsize=6)

    # ── right: energy ────────────────────────────────────────────────────────
    ax_en.set_xlim(0, t_axis[-1])
    ax_en.set_ylim(1 - 4e-14, 1 + 4e-14)
    ax_en.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax_en.set_xlabel("time $t$ (s)")
    ax_en.set_ylabel(r"$E_H(t)/E_H(0)$")
    ax_en.set_title("Discrete energy conservation")
    ax_en.axhline(1.0, color="0.6", lw=0.6, ls=":")
    line_en, = ax_en.plot([], [], lw=0.9, **{k: v for k, v in pubstyle.BW[0].items()
                                               if k != "marker"})

    fig.tight_layout()

    def update(i):
        tcf.set_array(snaps[i])
        time_txt.set_text(f"$t$ = {times[i]*1e6:.1f} µs")
        # energy up to this frame
        idx = i * SAVE_EVERY
        line_en.set_data(t_axis[:idx], E_hist[:idx])
        return tcf, time_txt, line_en

    ani = animation.FuncAnimation(fig, update, frames=n_frames,
                                  interval=50, blit=True)
    out = os.path.join(FIG, "anim_wave.gif")
    ani.save(out, writer=animation.PillowWriter(fps=20), dpi=60)
    plt.close(fig)
    print(f"  saved: {out}")


if __name__ == "__main__":
    make_animation()
