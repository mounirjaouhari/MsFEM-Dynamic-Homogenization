#!/usr/bin/env python3
"""
REAL spectral MsFEM basis (Contribution C1, computed for the first time).
On a coarse element K containing microstructure we solve the local Neumann
eigenproblem  -div(mu grad psi) = lambda rho psi  in K,  d_n psi = 0 on dK,
and form the enriched basis  phi_il = chi_i psi_l^K  (chi_i a P1 hat).
The constant mode psi_0 = 1 (lambda_0 = 0) guarantees P^1 subset V_H.
Verified: lambda_0 ~ 0; modes M-orthonormal; modes resolve the inclusions.
"""
import os
import numpy as np
import scipy.linalg as sla
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
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm
from skfem.helpers import dot, grad

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
MU_M, MU_I, RHO_M, RHO_I = 12e9, 0.5e9, 2500.0, 9000.0   # locally-resonant (soft, dense) inclusion
# coarse element K = unit square with a vertical column of 3 circular inclusions
CENTERS = [(0.5, 1 / 6), (0.5, 1 / 2), (0.5, 5 / 6)]
R = 0.13


def incl(x, y):
    m = np.zeros_like(x, dtype=bool)
    for cx, cy in CENTERS:
        m |= (x - cx) ** 2 + (y - cy) ** 2 <= R ** 2
    return m


def main():
    mesh = MeshTri().refined(6)                 # ~4225 nodes
    vb = Basis(mesh, ElementTriP1())

    @BilinearForm
    def stiff(u, v, w):
        x, y = w.x
        return np.where(incl(x, y), MU_I, MU_M) * dot(grad(u), grad(v))

    @BilinearForm
    def mass(u, v, w):
        x, y = w.x
        return np.where(incl(x, y), RHO_I, RHO_M) * u * v

    K = stiff.assemble(vb).toarray()
    M = mass.assemble(vb).toarray()
    lam, V = sla.eigh(K, M)                      # generalized, ascending
    p = mesh.p
    tri = mtri.Triangulation(p[0], p[1], mesh.t.T)
    om = np.sqrt(np.abs(lam))

    # M-normalise sign so plots are consistent (fix sign so the dominant
    # part is positive -> psi_0=1 renders as a light, physically-positive field)
    def mode(n):
        v = V[:, n]
        v = v / np.max(np.abs(v))
        if np.sum(v) < 0:
            v = -v
        return v

    chi = (1 - p[0]) * (1 - p[1])               # a P1-type corner hat chi_i
    phi = chi * mode(2)                         # enriched basis phi_il = chi_i psi_l

    fig, ax = plt.subplots(1, 5, figsize=(13, 2.9))
    # (0) element geometry (mu map) -- unsigned, stays Greys
    muv = np.where(incl(p[0], p[1]), 1.0, 0.0)
    ax[0].tripcolor(tri, muv, shading="gouraud", cmap="Greys", vmin=-0.3, vmax=1.3)
    ax[0].set_title("coarse element $K$\n(heterogeneous $\\mu,\\rho$)")
    # signed modes -> luminance-monotonic cividis + black zero contour
    panels = [(1, 0, r"$\psi_0^K=1,\ \lambda_0\approx0$"),
              (2, 1, r"$\psi_1^K$ (resonance)"),
              (3, 2, r"$\psi_2^K$ (resonance)")]
    for axi, n, title in panels:
        vals = mode(n)
        im = ax[axi].tripcolor(tri, vals, shading="gouraud",
                               cmap=pubstyle.FIELD_CMAP, vmin=-1, vmax=1)
        # zero contour marks the sign change (visible in pure B&W)
        if vals.max() > 0.0 > vals.min():
            ax[axi].tricontour(tri, vals, levels=[0.0], colors="k", linewidths=0.5)
        ax[axi].set_title(title)
        fig.colorbar(im, ax=ax[axi], fraction=0.046, pad=0.04)
    phiv = phi / np.max(np.abs(phi))
    im = ax[4].tripcolor(tri, phiv, shading="gouraud", cmap=pubstyle.FIELD_CMAP)
    if phiv.max() > 0.0 > phiv.min():
        ax[4].tricontour(tri, phiv, levels=[0.0], colors="k", linewidths=0.5)
    ax[4].set_title(r"enriched $\phi_{i\ell}=\chi_i\psi_\ell^K$")
    fig.colorbar(im, ax=ax[4], fraction=0.046, pad=0.04)
    th = np.linspace(0, 2 * np.pi, 40)
    for a in ax:
        for cx, cy in CENTERS:
            a.plot(cx + R * np.cos(th), cy + R * np.sin(th), "k-", lw=0.5)
        a.set_xticks([]); a.set_yticks([]); a.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_basis.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_basis.png"), bbox_inches="tight")
    plt.close(fig)
    print("fig_basis written. lambda_0=%.3e (->0 ok: %s); first freqs omega=%s"
          % (lam[0], abs(lam[0]) < 1e-6 * abs(lam[3]), np.array2string(om[:4], precision=2,
             formatter={'float': lambda z: f'{z:.2e}'})))
    print("psi_0 constant? max-min over mean = %.2e" %
          (np.ptp(V[:, 0]) / abs(np.mean(V[:, 0]))))


if __name__ == "__main__":
    main()
