#!/usr/bin/env python3
"""
Regenerate the last two Section-5 figures with REAL data:
  Fig 4 (sensitivity): homogenized modulus mu* vs contrast mu_i/mu_m
        (real cell-problem solves), shown inside the 2D Hashin-Shtrikman band.
  Fig 5 (wave profiles): a real computed wave-field snapshot U(x, y=1/2)
        from the verified transient solver.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skfem import MeshTri, Basis, ElementTriP1
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _d in ('solvers', 'utils'):
    _p = _os.path.join(_SRC, _d)
    if _p not in _sys.path: _sys.path.insert(0, _p)
del _sys, _os, _SRC, _d, _p
import cell_problem_2d as cpm

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
MU_M, VF = cpm.MU_M, cpm.VF


def mustar(contrast, refine=6):
    cpm.MU_I = contrast * MU_M
    m = MeshTri().refined(refine)
    vb = Basis(m, ElementTriP1())
    K = cpm.stiff_mu.assemble(vb)
    P, _, _ = cpm.periodic_prolongation(m)
    chi0 = cpm.solve_corrector(vb, K, P, cpm.load_dir(0).assemble(vb))
    ci = vb.interpolate(chi0)
    mu_bar = cpm.mu_avg.assemble(vb, chi=ci)
    return cpm.mu_times_dchi(0).assemble(vb, chi=ci) + mu_bar


def hs(contrast):
    mi = contrast * MU_M
    lo = MU_M + VF / (1.0 / (mi - MU_M) + (1 - VF) / (2 * MU_M))
    hi = mi + (1 - VF) / (1.0 / (MU_M - mi) + VF / (2 * mi))
    return min(lo, hi), max(lo, hi)


def sensitivity():
    cs = np.array([1.0001, 2, 3.12, 6.5, 10, 20, 40])
    ms = np.array([mustar(c) for c in cs]) / MU_M
    los, his = zip(*[hs(c) for c in cs])
    los = np.array(los) / MU_M; his = np.array(his) / MU_M
    fig, ax = plt.subplots(figsize=(4.3, 3.3))
    ax.fill_between(cs, los, his, color="0.85", label="Hashin--Shtrikman band")
    ax.loglog(cs, ms, "o-", color="#1f77b4", label=r"computed $\mu^*/\mu_m$")
    ax.set_xlabel(r"contrast $\mu_i/\mu_m$"); ax.set_ylabel(r"$\mu^*/\mu_m$")
    ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_sensitivity.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig 4 (sensitivity mu* vs contrast):")
    for c, m, lo, hi in zip(cs, ms, los, his):
        print(f"    contrast={c:6.2f}  mu*/mu_m={m:.3f}  HS=[{lo:.3f},{hi:.3f}]  in={lo-1e-3<=m<=hi+1e-3}")
    cpm.MU_I = 78.0e9   # restore


def wave():
    from skfem import BilinearForm
    from skfem.helpers import dot, grad
    MU, RHO = 1.747e10, 3825.0
    @BilinearForm
    def mass(u, v, w): return RHO * u * v
    @BilinearForm
    def stiff(u, v, w): return MU * dot(grad(u), grad(v))
    import scipy.sparse as sp, scipy.sparse.linalg as spla
    m = MeshTri().refined(6)
    vb = Basis(m, ElementTriP1())
    M = mass.assemble(vb); K = stiff.assemble(vb)
    free = vb.complement_dofs(vb.get_dofs())
    Mi = M[free][:, free]; Ki = K[free][:, free]
    ml = np.asarray(Mi.sum(axis=1)).ravel(); inv = 1.0 / ml
    lam = spla.eigsh(Ki, k=1, M=sp.diags(ml), which="LA", return_eigenvectors=False)[0]
    dt = 0.5 * 2.0 / np.sqrt(lam)
    x, y = m.p[0], m.p[1]
    # localized initial pulse near the left, zero velocity
    U0 = np.exp(-((x - 0.25)**2 + (y - 0.5)**2) / (2 * 0.06**2))
    Up = U0[free].copy()
    Uc = Up + 0.5 * dt**2 * (inv * (-(Ki @ Up)))
    nstep = 220
    for _ in range(nstep):
        Un = 2 * Uc - Up + dt**2 * (inv * (-(Ki @ Uc)))
        Up, Uc = Uc, Un
    U = np.zeros(K.shape[0]); U[free] = Uc
    # sample along y=1/2
    midline = np.isclose(y, 0.5, atol=1e-9)
    xs = x[midline]; us = U[midline]; o = np.argsort(xs)
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(xs[o], us[o], "-", color="#1f77b4", lw=1.4, label="computed field $U(x,0.5)$")
    ax.axvline(0.5, color="gray", ls=":", lw=0.9, label=r"interface $\Gamma$")
    ax.set_xlabel("position $x$"); ax.set_ylabel("displacement $U$")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_wave_profiles.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Fig 5 (wave snapshot): t={nstep*dt*1e3:.2f} ms, max|U|={np.max(np.abs(us)):.3f}")


if __name__ == "__main__":
    print("=" * 60)
    print("  Regenerating Fig 4 (sensitivity) and Fig 5 (wave)")
    print("=" * 60)
    sensitivity()
    wave()
    print("  done.")
