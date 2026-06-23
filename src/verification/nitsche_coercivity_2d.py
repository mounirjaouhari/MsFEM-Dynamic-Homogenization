#!/usr/bin/env python3
"""
GATE_VNV[Nitsche]: symmetric Nitsche / interior-penalty coercivity.

Assemble the SYMMETRIC interior-penalty (Nitsche) bilinear form for the 2D
Laplacian on a discontinuous P1 space,
    a_h(u,v) = sum_K (grad u, grad v)_K
             - sum_E <grad u . n>[v] - sum_E <grad v . n>[u]
             + (gamma0/h) sum_E [u][v]   (interior facets E + boundary),
and verify the discrete operator K_H(gamma0) is SYMMETRIC and becomes SPD
once gamma0 exceeds a threshold gamma_crit set by the inverse trace
inequality, while it is INDEFINITE below it. This is the coercivity claim
underlying Theorem (energy conservation) and the well-posedness.

K(gamma0) = A_bulk + A_cons + B_cons + gamma0 (A_pen + B_pen): linear in
gamma0, so we sweep cheaply and locate the SPD threshold.
"""
import numpy as np
from skfem import (MeshTri, Basis, ElementTriP1, ElementTriDG,
                   InteriorFacetBasis, FacetBasis, BilinearForm)
from skfem.helpers import grad, dot, jump


@BilinearForm
def bulk(u, v, w):
    return dot(grad(u), grad(v))


@BilinearForm
def cons_int(u, v, w):                 # -<grad u.n>[v] - <grad v.n>[u]
    n = w.n
    gu, gv = jump(w, dot(grad(u), n), dot(grad(v), n))
    ju, jv = jump(w, u, v)
    return -0.5 * (gu * jv + gv * ju)


@BilinearForm
def pen_int(u, v, w):                  # (1/h)[u][v]
    ju, jv = jump(w, u, v)
    return ju * jv / w.h


@BilinearForm
def cons_bnd(u, v, w):                 # boundary (Dirichlet via Nitsche)
    n = w.n
    return -dot(grad(u), n) * v - dot(grad(v), n) * u


@BilinearForm
def pen_bnd(u, v, w):
    return u * v / w.h


def symmetry_err(M):
    return abs((M - M.T)).max() / abs(M).max()


def main():
    print("=" * 64)
    print("  GATE_VNV[Nitsche] : symmetric interior-penalty coercivity")
    print("=" * 64)
    m = MeshTri().refined(3)
    e = ElementTriDG(ElementTriP1())
    ib = Basis(m, e)
    fb = InteriorFacetBasis(m, e)
    bb = FacetBasis(m, e)

    A = bulk.assemble(ib)
    Ac = cons_int.assemble(fb)
    Ap = pen_int.assemble(fb)
    Bc = cons_bnd.assemble(bb)
    Bp = pen_bnd.assemble(bb)

    K0 = (A + Ac + Bc)
    Kp = (Ap + Bp)
    n = K0.shape[0]
    print(f"  DG dofs = {n}")
    print(f"  symmetry: |K0-K0^T|/|K0| = {symmetry_err(K0.toarray()):.2e},  "
          f"|Kp-Kp^T|/|Kp| = {symmetry_err(Kp.toarray()):.2e}")

    gammas = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0]
    print("\n   gamma0     lambda_min        SPD?")
    lam_prev = None
    g_crit = None
    for g in gammas:
        K = (K0 + g * Kp).toarray()
        lam = np.linalg.eigvalsh(K)[0]
        spd = lam > 0
        print(f"   {g:6.1f}   {lam: .6e}    {spd}")
        if lam_prev is not None and lam_prev <= 0 < lam and g_crit is None:
            g_crit = (g_prev, g)
        lam_prev, g_prev = lam, g
    print("-" * 64)
    # refine threshold by bisection on the crossing bracket
    if g_crit:
        lo, hi = g_crit
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            lam = np.linalg.eigvalsh((K0 + mid * Kp).toarray())[0]
            if lam > 0:
                hi = mid
            else:
                lo = mid
        print(f"  gamma_crit ~ {hi:.3f}  (SPD for gamma0 >= gamma_crit)")
    ok = (symmetry_err(K0.toarray()) < 1e-10) and (g_crit is not None) \
        and (np.linalg.eigvalsh((K0 + 40.0 * Kp).toarray())[0] > 0)
    print(f"  GATE_VNV[Nitsche]: {'PASS' if ok else 'FAIL'} "
          f"(symmetric; indefinite below, SPD above gamma_crit)")


if __name__ == "__main__":
    main()
