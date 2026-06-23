#!/usr/bin/env python3
"""
Mesh figure (3 views), standard in the MsFEM/FEM literature:
 (a) coarse MsFEM mesh T_H with the mean-line interface Gamma (the online solve);
 (b) fine mesh resolving the row of circular fibres (the inclusion-resolving
     reference);
 (c) zoom on the enlarged interface of thickness e (faces Gamma^- , Gamma^+).
Shows the scale separation h << H and the actual discretisations used.
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
import matplotlib.tri as mtri
from skfem import MeshTri

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
W = 2.0
P = 0.25                      # transverse period of the row
R = 0.09                      # fibre radius
E = 0.20                      # enlarged-interface thickness e = O(h)
CENTERS_Y = np.arange(P / 2, 1.0, P)


def tri_of(nx, ny):
    m = MeshTri.init_tensor(np.linspace(-W, W, nx), np.linspace(0, 1, ny))
    return mtri.Triangulation(m.p[0], m.p[1], m.t.T)


def draw_fibres(ax, lw=0.8, fill=True):
    th = np.linspace(0, 2 * np.pi, 60)
    for cy in CENTERS_Y:
        if fill:
            ax.fill(R * np.cos(th), cy + R * np.sin(th), color="0.7", alpha=0.6, zorder=3)
        ax.plot(R * np.cos(th), cy + R * np.sin(th), "k-", lw=lw, zorder=4)


def main():
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
    # (a) coarse MsFEM mesh + interface Gamma
    ax[0].triplot(tri_of(9, 5), color="#1f77b4", lw=0.8)
    ax[0].axvline(0, color="#c0392b", lw=2.4, zorder=5)
    ax[0].text(0.14, 1.03, r"$\Gamma$", color="#c0392b", fontsize=12)
    ax[0].annotate("", xy=(-2.0, -0.09), xytext=(-1.5, -0.09),
                   arrowprops=dict(arrowstyle="<->", color="k"))
    ax[0].text(-1.82, -0.2, "$H$", fontsize=11)
    ax[0].set_title("(a) coarse MsFEM mesh $\\mathcal{T}_H$\n(online solve, bulk = matrix)")
    # (b) fine resolved mesh + fibres
    ax[1].triplot(tri_of(121, 31), color="0.6", lw=0.25)
    draw_fibres(ax[1], lw=0.7)
    ax[1].annotate("", xy=(0.55, CENTERS_Y[1]), xytext=(0.55, CENTERS_Y[2]),
                   arrowprops=dict(arrowstyle="<->", color="#1a8a3a", lw=1.6))
    ax[1].text(0.62, 0.5 * (CENTERS_Y[1] + CENTERS_Y[2]), "$h$", color="#1a8a3a", fontsize=12)
    ax[1].set_title("(b) fine mesh resolving the row\n(inclusion-resolving reference)")
    # (c) zoom on enlarged interface
    ax[2].triplot(tri_of(121, 31), color="0.6", lw=0.5)
    draw_fibres(ax[2], lw=1.0)
    ax[2].axvspan(-E / 2, E / 2, facecolor="#f9d29b", edgecolor="none",
                  hatch="////", alpha=0.5, zorder=2)
    for s, lab in ((-E / 2, r"$\Gamma^-$"), (E / 2, r"$\Gamma^+$")):
        ax[2].axvline(s, color="#d68910", ls="--", lw=1.3, zorder=5)
    ax[2].axvline(0, color="#c0392b", lw=1.8, zorder=5)
    ax[2].annotate("", xy=(-E / 2, 0.50), xytext=(E / 2, 0.50),
                   arrowprops=dict(arrowstyle="<->", color="#b9770e", lw=1.4))
    ax[2].text(0.065, 0.52, "$e$", color="#b9770e", fontsize=11, ha="center")
    ax[2].text(-E / 2 - 0.02, 0.70, r"$\Gamma^-$", color="#d68910", fontsize=11, ha="right")
    ax[2].text(E / 2 + 0.02, 0.70, r"$\Gamma^+$", color="#d68910", fontsize=11, ha="left")
    ax[2].set_xlim(-0.42, 0.42); ax[2].set_ylim(0.25, 0.75)
    ax[2].set_title("(c) enlarged interface, thickness $e=\\mathcal{O}(h)$\n(faces $\\Gamma^\\pm$ at $\\pm e/2$)")
    for a in ax:
        a.set_aspect("equal"); a.set_xlabel("$x$"); a.set_ylabel("$y$")
    ax[0].set_ylim(-0.2, 1.12)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_mesh.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_mesh.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig_mesh written (coarse + fine + zoom); period h=%.2f, R=%.2f, e=%.2f" % (P, R, E))


if __name__ == "__main__":
    main()
