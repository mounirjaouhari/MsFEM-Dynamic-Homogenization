#!/usr/bin/env python3
"""Generate the Model-A Section-5 figures from the verified solvers."""
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
import transmission_resolved as ref
import transmission_compare as tc
import transmission_oblique as to
import resonant_mass_2d as rm

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)


def fig_gateref():
    W = 6.0
    etas = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    eb, ej, Be, Se = [], [], [], []
    for eta in etas:
        tr, rr, _, _ = ref.solve(W, "circle", eta)
        tj, _ = tc.tau_jump(eta)
        eb.append(abs(tr - 1.0)); ej.append(abs(tr - tj))
        be, se = tc.invert(tr, rr, eta)
        Be.append(be); Se.append(se)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.6, 5.6))
    ax1.loglog(etas, eb, label=r"no model: $|\tau_{\rm res}-1|$", **pubstyle.BW[0])
    ax1.loglog(etas, ej, label=r"jump model: $|\tau_{\rm res}-\tau_{\rm jump}|$",
               **pubstyle.BW[1])
    ax1.loglog(etas, 0.1 * etas, ls=":", color="0.0", lw=0.8,
               label=r"$\mathcal{O}(\eta)$")
    ax1.set_xlabel(r"$\eta=kh$"); ax1.set_ylabel("transmission error")
    ax1.legend(); ax1.grid(True, which="both", alpha=0.3)
    ax1.set_title("jump model reproduces the\nresolved row")
    ax2.semilogx(etas, Be, label=r"$\mathcal{B}_{\rm eff}$ (inverted)", **pubstyle.BW[0])
    ax2.semilogx(etas, Se, label=r"$\mathcal{S}_{\rm eff}$ (inverted)", **pubstyle.BW[1])
    ax2.axhline(tc.B, ls=":", color="0.0", lw=0.8, label=r"cell $\mathcal{B}=%.2f$" % tc.B)
    ax2.axhline(tc.S, ls=(0, (1, 1)), color="0.45", lw=0.8,
                label=r"cell $\mathcal{S}=%.2f$" % tc.S)
    ax2.set_xlabel(r"$\eta=kh$"); ax2.set_ylabel("effective parameter")
    ax2.legend(); ax2.grid(True, alpha=0.3)
    ax2.set_title("two-route parameter\nagreement")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_gateref.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_gateref.png"), bbox_inches="tight")
    plt.close(fig)
    print("  fig_gateref.pdf/png  (last eta: B_eff=%.3f S_eff=%.3f vs cell %.3f/%.3f)"
          % (Be[-1], Se[-1], tc.B, tc.S))


def fig_bandgap():
    omega_n, a_n, M_incl, _ = rm.inclusion_modes()
    w1 = omega_n[0]
    wg = np.linspace(0.05 * w1, 1.6 * w1, 4000)
    Me = rm.effective_mass(wg, omega_n, a_n, M_incl) / M_incl
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.axhline(0, color="k", lw=0.6)
    neg = Me < 0
    ax.fill_between(wg / w1, -6, 6, where=neg, facecolor="#f4cccc", edgecolor="0.4",
                    hatch="////", linewidth=0.5,
                    label="band gap ($\\mathcal{S}_{\\rm eff}<0$)")
    ax.plot(wg / w1, np.clip(Me, -6, 6), color="#0072B2", ls="-", lw=1.6,
            label=r"$\mathcal{S}_{\rm eff}(\omega)/\mathcal{S}_0$")
    ax.axvline(1.0, ls="--", color="#D55E00", lw=1.0, label=r"resonance $\omega_1$")
    ax.set_xlabel(r"$\omega/\omega_1$"); ax.set_ylabel(r"$\mathcal{S}_{\rm eff}(\omega)/\mathcal{S}_0$")
    ax.set_ylim(-6, 6); ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_bandgap.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_bandgap.png"), bbox_inches="tight")
    plt.close(fig)
    gap = wg[neg]
    print("  fig_bandgap.pdf/png  (gap %.2f-%.2f x omega_1; f1=%.1f kHz)"
          % (gap.min() / w1, gap.max() / w1, w1 / 2 / np.pi / 1e3))


def fig_field2d():
    """Genuinely-2D worked example: resolved oblique scattering through the row."""
    W, eta, theta = 4.0, 0.45, np.deg2rad(45.0)
    tau, refl, U, mesh, (kx, ky) = to.solve_resolved_oblique(
        W, "circle", eta, theta, n_in=181, n_out=30, ny=61, return_field=True)
    p = mesh.p
    nrep = 4                       # stack periods via Bloch phase e^{i n ky}
    xs, ys, vs, tris = [], [], [], []
    base = mesh.t.T
    for n in range(nrep):
        xs.append(p[0]); ys.append(p[1] + n)
        vs.append(np.real(np.exp(1j * n * ky) * U))
        tris.append(base + n * p.shape[1])
    X = np.concatenate(xs); Y = np.concatenate(ys)
    V = np.concatenate(vs); T = np.concatenate(tris, axis=0)
    tri = mtri.Triangulation(X, Y, T)
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    lim = 0.9 * np.max(np.abs(V))
    tpc = ax.tripcolor(tri, V, shading="gouraud", cmap=pubstyle.FIELD_CMAP,
                       vmin=-lim, vmax=lim)
    ax.tricontour(tri, V, levels=[0], colors="k", linewidths=0.5)  # sign change in B&W
    th = np.linspace(0, 2 * np.pi, 80)
    for n in range(nrep):                      # inclusion outlines at x=0
        ax.plot(ref.R * np.cos(th), 0.5 + n + ref.R * np.sin(th), "k-", lw=0.6)
    ax.axvline(0, color="k", ls=":", lw=0.7)
    ax.set_xlim(-W, W); ax.set_ylim(0, nrep)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$ (stacked periods)")
    ax.set_title(r"resolved 2D field $\mathrm{Re}\,U(x,y)$, $45^\circ$ incidence")
    fig.colorbar(tpc, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_field2d.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_field2d.png"), bbox_inches="tight")
    plt.close(fig)
    print("  fig_field2d.pdf/png  (eta=%.2f, 45deg, |tau|=%.4f) -- genuine 2D wavefield"
          % (eta, abs(tau)))


if __name__ == "__main__":
    print("Generating Model-A figures ->", os.path.normpath(FIG))
    fig_gateref()
    fig_bandgap()
    fig_field2d()
    print("done.")
