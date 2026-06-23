#!/usr/bin/env python3
"""
Six-window solution panel; BOTH the resolved row and the homogenized model are
REAL finite-element solves (no drawn data).  The static corrector (unit normal
shear gradient) gives a real, large displacement jump [U]=h*B, ideal to see.
 Top:  (a) resolved row U-x (continuous, varies through the layer),
       (b) homogenized FE solve U-x (matrix + Nitsche jump -> sharp step),
       (c) resolved dynamic field at 45 deg (genuinely 2D).
 Bottom:(d) cut U-x across Gamma (resolved vs homogenized FE solve),
       (e) zoom on the jump [U], (f) jump B vs phase contrast.
Fields centred (subtract mean) since the pure-Neumann corrector is defined up to
a constant.  Sources: interface_params_2d, homog_fe, transmission_oblique.
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
import homog_fe as hf
import transmission_oblique as to

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
L = 4.0


def tcol(ax, mesh, vals, vlim, ystack=1, ky=0.0, cmap=pubstyle.FIELD_CMAP, zero=True):
    """Plot a (possibly stacked) field with cividis; overlay a black zero contour
    for signed fields so the sign change is visible in pure grayscale."""
    p = mesh.p
    xs, ys, vv, tt = [], [], [], []
    base = mesh.t.T
    for n in range(ystack):
        xs.append(p[0]); ys.append(p[1] + n)
        vv.append(np.real(np.exp(1j * n * ky) * vals)); tt.append(base + n * p.shape[1])
    X = np.concatenate(xs); Y = np.concatenate(ys); V = np.concatenate(vv)
    T = np.concatenate(tt, axis=0)
    tri = mtri.Triangulation(X, Y, T)
    tc = ax.tripcolor(tri, V, shading="gouraud", cmap=cmap, vmin=-vlim, vmax=vlim)
    if zero:
        ax.tricontour(tri, V, levels=[0], colors="k", linewidths=0.5)
    return tc


def yavg(mesh, w):
    p = mesh.p.T
    xs = np.unique(np.round(p[:, 0], 9))
    return xs, np.array([w[np.abs(p[:, 0] - x) < 1e-9].mean() for x in xs])


def static_corrector(contrast=None):
    if contrast is not None:
        ip.MU_I = contrast * ip.MU_M
    mesh = ip.build_strip(L, nx=min(8 * 60 + 1, int(2 * L * 60) + 1), ny=61)
    mu = ip.make_mu("circle")
    vb, K = ip.assemble(mesh, mu)
    u, B, C1 = ip.solve_normal(mesh, vb, K, mu, L)
    return mesh, u - mesh.p[0], B


def main():
    print("Six-window jump panel (both panels real FE solves) ->", os.path.normpath(FIG))
    ip.MU_I = 78.0e9
    meshR, wR, B = static_corrector()
    wR = wR - wR.mean()                                   # centre (Neumann corrector)
    hr = hf.solve_homog_static(B, L=L, nx=120, ny=61)     # REAL homogenized FE solve
    wHL = hr["UL"] - hr["mL"].p[0]; wHR = hr["UR"] - hr["mR"].p[0]
    off = 0.5 * (wHL.mean() + wHR.mean()); wHL -= off; wHR -= off
    tau_o, _, Uo, mesho, (kxo, kyo) = to.solve_resolved_oblique(
        L, "circle", 0.45, np.deg2rad(45), n_in=181, n_out=30, ny=61, return_field=True)

    vlim = 0.85 * max(np.abs(wR).max(), np.abs(wHL).max(), np.abs(wHR).max())
    fig, ax = plt.subplots(2, 3, figsize=(11, 6))

    # (a) resolved row, real FE
    th = np.linspace(0, 2 * np.pi, 60)
    tc = tcol(ax[0, 0], meshR, wR, vlim, ystack=3)
    for n in range(3):
        ax[0, 0].plot(ip.R * np.cos(th), 0.5 + n + ip.R * np.sin(th), "k-", lw=0.5)
    ax[0, 0].axvline(0, color="k", ls=":", lw=0.6)
    ax[0, 0].set_title("(a) resolved row $U-x$ (FE)\ncontinuous through the layer")
    fig.colorbar(tc, ax=ax[0, 0], shrink=0.8)
    # (b) homogenized FE solve (two half-strips)
    tcol(ax[0, 1], hr["mL"], wHL, vlim, ystack=3)
    tcb = tcol(ax[0, 1], hr["mR"], wHR, vlim, ystack=3)
    ax[0, 1].axvline(0, color="k", lw=1.2)
    ax[0, 1].set_title("(b) homogenized FE solve $U-x$\nmatrix + Ventcel jump $[U]=h\\mathcal{B}$")
    fig.colorbar(tcb, ax=ax[0, 1], shrink=0.8)
    # (c) resolved dynamic 45 deg
    tcd = tcol(ax[0, 2], mesho, Uo, 0.85 * np.abs(np.real(Uo)).max(), ystack=3, ky=kyo)
    ax[0, 2].axvline(0, color="k", ls=":", lw=0.6)
    ax[0, 2].set_title("(c) resolved dynamic $45^\\circ$ (FE)\ngenuinely 2D, $\\partial_yU\\neq0$")
    fig.colorbar(tcd, ax=ax[0, 2], shrink=0.8)
    for a in ax[0]:
        a.set_xlim(-L, L); a.set_ylim(0, 3); a.set_xlabel("$x$"); a.set_ylabel("$y$")

    # (d) cut U-x across Gamma: resolved vs homogenized FE
    xs, wx = yavg(meshR, wR)
    xL, wlx = yavg(hr["mL"], wHL); xR, wrx = yavg(hr["mR"], wHR)
    me0 = max(1, len(xs) // 8); me1 = max(1, len(xL) // 8)
    ax[1, 0].plot(xs, wx, **pubstyle.BW[0], markevery=me0, label="resolved (FE)")
    ax[1, 0].plot(xL, wlx, **pubstyle.BW[1], markevery=me1, label="homogenized (FE)")
    ax[1, 0].plot(xR, wrx, **pubstyle.BW[1], markevery=me1)
    ax[1, 0].axvline(0, color="k", ls=":", lw=0.6)
    ax[1, 0].set_title("(d) cut $U-x$ across $\\Gamma$")
    ax[1, 0].set_xlabel("$x$"); ax[1, 0].legend(); ax[1, 0].grid(alpha=0.3)
    # (e) zoom on the jump
    m = np.abs(xs) < 1.2
    mL = np.abs(xL) < 1.2; mR = np.abs(xR) < 1.2
    meZ = max(1, int(m.sum()) // 8); meZL = max(1, int(mL.sum()) // 8)
    ax[1, 1].plot(xs[m], wx[m], **pubstyle.BW[0], markevery=meZ, label="resolved")
    ax[1, 1].plot(xL[mL], wlx[mL], **pubstyle.BW[1], markevery=meZL, label="homogenized")
    ax[1, 1].plot(xR[mR], wrx[mR], **pubstyle.BW[1], markevery=meZL)
    jH = float(np.interp(0, xR, wrx) - np.interp(0, xL[::-1], wlx[::-1]))
    ax[1, 1].annotate("", xy=(0, np.interp(0, xR, wrx)), xytext=(0, np.interp(0, xL[::-1], wlx[::-1])),
                      arrowprops=dict(arrowstyle="<->", color="k", lw=1.8))
    ax[1, 1].text(0.07, 0, r"$[U]=h\mathcal{B}$", color="k", fontsize=10)
    ax[1, 1].axvline(0, color="k", ls=":", lw=0.6)
    ax[1, 1].set_title("(e) zoom: displacement jump $[U]$")
    ax[1, 1].set_xlabel("$x$"); ax[1, 1].legend(); ax[1, 1].grid(alpha=0.3)
    # (f) B vs contrast
    contrasts = np.array([1.5, 2, 3, 6.5, 12, 25]); Bs = []
    for c in contrasts:
        _, _, Bc = static_corrector(c); Bs.append(Bc)
    ip.MU_I = 78.0e9
    ax[1, 2].semilogx(contrasts, Bs, **pubstyle.BW[0])
    ax[1, 2].axhline(0, color="k", ls=":", lw=0.8)
    ax[1, 2].set_title("(f) jump law: $\\mathcal{B}$ vs contrast")
    ax[1, 2].set_xlabel(r"$\mu_i/\mu_m$"); ax[1, 2].set_ylabel(r"$\mathcal{B}$"); ax[1, 2].grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_panel6.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_panel6.png"), bbox_inches="tight")
    plt.close(fig)
    print("  B(circle)=%.3f  homog-FE jump=%.3f  zoom-jump=%.3f  oblique|tau|=%.4f" %
          (B, hr["jump"], jH, abs(tau_o)))
    print("  centred resolved w[%.3f,%.3f]  homog w[%.3f,%.3f]" %
          (wR.min(), wR.max(), min(wHL.min(), wHR.min()), max(wHL.max(), wHR.max())))


if __name__ == "__main__":
    main()
