#!/usr/bin/env python3
"""
Real strip-cell corrector fields that DEFINE the interface parameters (Model A).
Left: normal corrector Q^(1) = U - x (far-field jump = B, the displacement-jump
compliance).  Right: tangential corrector Q^(2) (its weighted tangential flux
int mu d_y2 Q^(2) = C; for a centred circle the far-field jump B2 ~ 0 by symmetry).
Both are verified solves of the static strip cell problem (interface_params_2d).
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
import interface_params_2d as ip

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
L = 4.0


def tcol(ax, mesh, vals, vlim, ystack=3, cmap=pubstyle.FIELD_CMAP):
    p = mesh.p
    xs, ys, vv, tt = [], [], [], []
    base = mesh.t.T
    for n in range(ystack):
        xs.append(p[0]); ys.append(p[1] + n); vv.append(vals); tt.append(base + n * p.shape[1])
    X = np.concatenate(xs); Y = np.concatenate(ys); V = np.concatenate(vv); T = np.concatenate(tt, axis=0)
    tri = mtri.Triangulation(X, Y, T)
    tpc = ax.tripcolor(tri, V, shading="gouraud", cmap=cmap, vmin=-vlim, vmax=vlim)
    # signed field: black zero contour so the sign change reads in B&W
    ax.tricontour(tri, V, levels=[0], colors="k", linewidths=0.5)
    th = np.linspace(0, 2 * np.pi, 60)
    for n in range(ystack):
        ax.plot(ip.R * np.cos(th), 0.5 + n + ip.R * np.sin(th), "k-", lw=0.6)
    ax.axvline(0, color="k", ls=":", lw=0.6)
    ax.set_xlim(-2.6, 2.6); ax.set_ylim(0, ystack); ax.set_xlabel("$y_1$"); ax.set_ylabel("$y_2$")
    return tpc


def main():
    ip.MU_I = 78.0e9
    mesh = ip.build_strip(L, nx=min(8 * 60 + 1, int(2 * L * 60) + 1), ny=61)
    mu = ip.make_mu("circle")
    vb, K = ip.assemble(mesh, mu)
    u, B, C1 = ip.solve_normal(mesh, vb, K, mu, L)
    W, B2, C = ip.solve_tangential(mesh, vb, K, mu, L)
    Q1 = u - mesh.p[0]; Q1 = Q1 - Q1.mean()
    Q2 = W - W.mean()
    v1 = 0.85 * np.abs(Q1).max(); v2 = 0.85 * np.abs(Q2).max()

    fig, ax = plt.subplots(2, 1, figsize=(5.2, 6.0))
    tc1 = tcol(ax[0], mesh, Q1, v1)
    ax[0].set_title(r"(a) normal corrector $Q^{(1)}=U-y_1$" + "\n"
                    + r"far-field jump $=\mathcal{B}=%.2f$" % B, fontsize=9)
    fig.colorbar(tc1, ax=ax[0], shrink=0.85)
    # annotate the plateau jump (black so it reads in B&W)
    ax[0].annotate("", xy=(2.3, 0.5 + 1 - 0.0), xytext=(2.3, 0.5 + 1 + B),
                   arrowprops=dict(arrowstyle="<->", color="k", lw=1.4))
    tc2 = tcol(ax[1], mesh, Q2, v2)
    ax[1].set_title(r"(b) tangential corrector $Q^{(2)}$" + "\n"
                    + r"$\mathcal{C}/\mu_m=%.2f$, $\mathcal{B}_2\approx0$" % (C / ip.MU_M), fontsize=9)
    fig.colorbar(tc2, ax=ax[1], shrink=0.85)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_cellcorr.pdf"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_cellcorr.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("fig_cellcorr written; B=%.3f  B2=%.2e  C/mu_m=%.3f  C1/mu_m=%.2e"
          % (B, B2, C / ip.MU_M, C1 / ip.MU_M))


if __name__ == "__main__":
    main()
