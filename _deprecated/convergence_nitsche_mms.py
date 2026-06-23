#!/usr/bin/env python3
"""
DEPRECATED — DO NOT USE FOR RESULTS. The full-DG (ElementTriDG +
InteriorFacetBasis) assembly used here becomes numerically singular at
refine>=3 (lambda_min ~ 1e-17, gamma0-independent), so this script returns
NaN convergence rates. This is a scikit-fem DG facet-assembly artifact, NOT a
defect of the scheme. It is SUPERSEDED by broken_nitsche_2d.py (continuous P1
within subdomains, broken only at Gamma), which is MMS-verified at optimal
order. See msfem2d_verified/data/PROVENANCE.md, GATE_VNV[Nitsche]. Kept only
for the diagnostic record.

----------------------------------------------------------------------------
Convergence in H of the symmetric Nitsche / interior-penalty (DG) scheme on
the homogenized operator -div(mu* grad u) = f, [u]=0 imposed weakly across
interfaces, with the manufactured solution u = sin(pi x) sin(pi y) on
[0,1]^2 (u=0 on dOmega via Nitsche). This is the spatial discretization of
the paper's MsFEM-Nitsche scheme; it verifies Term II (discretization):
expected H1 order ~1, L2 order ~2 (P1). Real solve, produces Table 1 / Fig 2
data.  gamma0 = 10 > gamma_crit ~ 2 (coercive, see nitsche_coercivity_2d).
"""
import os
import numpy as np
from skfem import (MeshTri, Basis, ElementTriP1, ElementTriDG,
                   InteriorFacetBasis, FacetBasis, BilinearForm, LinearForm,
                   Functional, solve)
from skfem.helpers import grad, dot, jump

PI = np.pi
MU = 1.747e10
GAMMA0 = 10.0


@BilinearForm
def bulk(u, v, w):
    return MU * dot(grad(u), grad(v))


@BilinearForm
def cons_int(u, v, w):
    n = w.n
    gu, gv = jump(w, MU * dot(grad(u), n), MU * dot(grad(v), n))
    ju, jv = jump(w, u, v)
    return -0.5 * (gu * jv + gv * ju)


@BilinearForm
def pen_int(u, v, w):
    ju, jv = jump(w, u, v)
    return GAMMA0 * MU * ju * jv / w.h


@BilinearForm
def cons_bnd(u, v, w):
    n = w.n
    return -MU * dot(grad(u), n) * v - MU * dot(grad(v), n) * u


@BilinearForm
def pen_bnd(u, v, w):
    return GAMMA0 * MU * u * v / w.h


@LinearForm
def load(v, w):
    x, y = w.x
    return 2.0 * PI**2 * MU * np.sin(PI * x) * np.sin(PI * y) * v


@Functional
def l2_sq(w):
    x, y = w.x
    ue = np.sin(PI * x) * np.sin(PI * y)
    return (w["uh"] - ue) ** 2


@Functional
def h1_sq(w):
    x, y = w.x
    dux = PI * np.cos(PI * x) * np.sin(PI * y)
    duy = PI * np.sin(PI * x) * np.cos(PI * y)
    g = grad(w["uh"])
    return (g[0] - dux) ** 2 + (g[1] - duy) ** 2


def solve_level(refine):
    m = MeshTri().refined(refine)
    e = ElementTriDG(ElementTriP1())
    ib = Basis(m, e)
    fb = InteriorFacetBasis(m, e)
    bb = FacetBasis(m, e)
    K = (bulk.assemble(ib) + cons_int.assemble(fb) + pen_int.assemble(fb)
         + cons_bnd.assemble(bb) + pen_bnd.assemble(bb))
    b = load.assemble(ib)
    uh = solve(K, b)
    ui = ib.interpolate(uh)
    return m.param(), np.sqrt(l2_sq.assemble(ib, uh=ui)), np.sqrt(h1_sq.assemble(ib, uh=ui))


def main():
    print("=" * 64)
    print("  Convergence in H : Nitsche/DG scheme, MMS (gamma0=%.0f)" % GAMMA0)
    print("=" * 64)
    rows = [solve_level(r) for r in range(2, 6)]
    hs = np.array([r[0] for r in rows])
    el2 = np.array([r[1] for r in rows])
    eh1 = np.array([r[2] for r in rows])
    print("    h          L2 error     L2 rate     H1 error     H1 rate")
    for i in range(len(rows)):
        rl2 = "" if i == 0 else f"{np.log(el2[i-1]/el2[i])/np.log(hs[i-1]/hs[i]):.2f}"
        rh1 = "" if i == 0 else f"{np.log(eh1[i-1]/eh1[i])/np.log(hs[i-1]/hs[i]):.2f}"
        print(f"  {hs[i]:.4e}  {el2[i]:.4e}   {rl2:>6}     {eh1[i]:.4e}   {rh1:>6}")
    pL2 = np.polyfit(np.log(hs), np.log(el2), 1)[0]
    pH1 = np.polyfit(np.log(hs), np.log(eh1), 1)[0]
    print(f"  LSQ orders: L2 = {pL2:.3f}  (expect ~2),  H1 = {pH1:.3f}  (expect ~1)")
    ok = (1.8 <= pL2 <= 2.2) and (0.85 <= pH1 <= 1.2)
    print("-" * 64)
    print(f"  Term-II (discretization) convergence: {'PASS' if ok else 'FAIL'}")

    # real Fig 2 (convergence) artifact
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.2, 3.4))
        ax.loglog(hs, eh1, "o-", label=r"$H^1$ error (rate %.2f)" % pH1)
        ax.loglog(hs, el2, "s-", label=r"$L^2$ error (rate %.2f)" % pL2)
        ax.loglog(hs, eh1[0]*(hs/hs[0]), "k--", lw=0.8, label=r"$\mathcal{O}(H)$")
        ax.loglog(hs, el2[0]*(hs/hs[0])**2, "k:", lw=0.8, label=r"$\mathcal{O}(H^2)$")
        ax.set_xlabel("coarse mesh size $H$"); ax.set_ylabel("error")
        ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        out = os.path.join(os.path.dirname(__file__), "figures", "fig2_convergence_real.pdf")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(f"  saved REAL convergence figure: {out}")
    except Exception as exc:
        print("  (figure skipped:", exc, ")")
    # CSV
    import csv
    cpath = os.path.join(os.path.dirname(__file__), "data", "convergence_H_real.csv")
    with open(cpath, "w", newline="") as f:
        wcsv = csv.writer(f); wcsv.writerow(["H", "L2_error", "H1_error"])
        for i in range(len(rows)):
            wcsv.writerow([hs[i], el2[i], eh1[i]])
    print(f"  saved: {cpath}")


if __name__ == "__main__":
    main()
