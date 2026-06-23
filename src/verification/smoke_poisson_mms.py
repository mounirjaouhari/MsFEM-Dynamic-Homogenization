#!/usr/bin/env python3
"""
GATE_SMOKE (workflow S2.0) -- environment smoke test.

Solve the trivial 2D Poisson problem  -Laplace(u) = f  on the unit square,
u = 0 on the boundary, with the manufactured solution
    u_exact(x,y) = sin(pi x) sin(pi y),  f = 2 pi^2 sin(pi x) sin(pi y),
and verify the FEM stack (scikit-fem, P1 triangles) converges at the
expected order: O(h^2) in L2, O(h) in H1. Also check determinism.

This is a REAL computation (no hardcoded results). Pass criterion:
observed L2 order in [1.85, 2.15], H1 order in [0.85, 1.15], deterministic.
"""
import numpy as np
from skfem import (MeshTri, Basis, ElementTriP1, BilinearForm, LinearForm,
                   Functional, solve, condense)
from skfem.helpers import dot, grad

PI = np.pi


@BilinearForm
def stiffness(u, v, w):
    return dot(grad(u), grad(v))


@LinearForm
def load(v, w):
    x, y = w.x
    return 2.0 * PI**2 * np.sin(PI * x) * np.sin(PI * y) * v


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
    vh = Basis(m, ElementTriP1())
    A = stiffness.assemble(vh)
    b = load.assemble(vh)
    D = vh.get_dofs()                       # boundary DOFs (Dirichlet u=0)
    uh = solve(*condense(A, b, D=D))
    ui = vh.interpolate(uh)
    e_l2 = np.sqrt(l2_sq.assemble(vh, uh=ui))
    e_h1 = np.sqrt(h1_sq.assemble(vh, uh=ui))
    h = m.param()                           # representative mesh size
    return h, e_l2, e_h1, vh.N


def run():
    rows = [solve_level(r) for r in range(2, 7)]
    hs = np.array([r[0] for r in rows])
    el2 = np.array([r[1] for r in rows])
    eh1 = np.array([r[2] for r in rows])
    ndof = [r[3] for r in rows]
    # least-squares slopes of log(err) vs log(h)
    p_l2 = np.polyfit(np.log(hs), np.log(el2), 1)[0]
    p_h1 = np.polyfit(np.log(hs), np.log(eh1), 1)[0]
    print("  level  ndof      h         L2err        H1err")
    for r in range(len(rows)):
        print(f"  {r+2:>5d} {ndof[r]:>6d}  {hs[r]:.4e}  {el2[r]:.4e}  {eh1[r]:.4e}")
    print(f"  observed orders:  L2 = {p_l2:.3f}   H1 = {p_h1:.3f}")
    return p_l2, p_h1, el2


if __name__ == "__main__":
    print("=" * 64)
    print("  GATE_SMOKE : 2D Poisson MMS (scikit-fem, P1)")
    print("=" * 64)
    p_l2, p_h1, el2 = run()
    # determinism: identical results on a repeat run
    el2_b = run()[2]
    deterministic = bool(np.array_equal(el2, el2_b))
    ok = (1.85 <= p_l2 <= 2.15) and (0.85 <= p_h1 <= 1.15) and deterministic
    print("-" * 64)
    print(f"  deterministic re-run: {deterministic}")
    print(f"  GATE_SMOKE: {'PASS' if ok else 'FAIL'}  "
          f"(L2~2? {1.85<=p_l2<=2.15}; H1~1? {0.85<=p_h1<=1.15})")
