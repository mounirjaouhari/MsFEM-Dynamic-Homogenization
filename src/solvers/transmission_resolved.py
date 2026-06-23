#!/usr/bin/env python3
"""
MODEL A validation, part 1/2: the INCLUSION-RESOLVING REFERENCE (GATE_REF).
Time-harmonic antiplane shear transmission across a periodic row of inclusions,
normal incidence.  This is the fine-scale truth the homogenized-jump model
(transmission_homogenized.py) must reproduce to O(eta^2).

Geometry: x in (-W, W), y in (0, 1) (one period, 1-periodic in y), one circular
inclusion (mu_i, rho_i, radius R) centred at (0, 1/2); matrix (mu_m, rho_m)
elsewhere.  Period h = 1, so eta = k_m * h = k_m.

PDE:  div(mu grad U) + rho omega^2 U = 0  (time-harmonic, e^{-i omega t}).
Incident plane wave from the left  U_inc = e^{i k_m (x+W)}.
First-order ABC (EXACT for the propagating 0th order at normal incidence):
  right (x=+W, outgoing):  d_x U = i k_m U
  left  (x=-W):            d_x U = -i k_m U + 2 i k_m U_inc
Weak form:  a(U,V) = int mu gradU.gradV - omega^2 int rho U V
                     - i k_m mu_m ( int_{left}UV + int_{right}UV )
            b(V)   = -2 i k_m mu_m int_{left} U_inc V .
Observables:  transmission  tau = <U>_{x=+W} e^{-i k_m 2W};  reflection rho = <U>_{x=-W} - 1.

VERIFICATION: matrix-only (no inclusion) must give |tau|=1, |rho|=0 EXACTLY
(the plane wave solves the discrete system since the ABC is exact for it).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import (MeshTri, Basis, ElementTriP1, FacetBasis,
                   BilinearForm, LinearForm)
from skfem.helpers import dot, grad
from interface_params_2d import y_periodic_prolongation, MU_M, MU_I, RHO_M, RHO_I, R


def graded_x(W, n_in, n_out, x_in=2.0):
    """Dense nodes in (-x_in, x_in) around the row, coarser out to +-W."""
    inner = np.linspace(-x_in, x_in, n_in)
    if W > x_in:
        rt = x_in + (W - x_in) * (np.linspace(0, 1, n_out) ** 1.0)
        x = np.unique(np.concatenate([-rt[::-1], inner, rt]))
    else:
        x = inner
    return x


def build(W, n_in=161, n_out=40, ny=49):
    x = graded_x(W, n_in, n_out)
    y = np.linspace(0.0, 1.0, ny)
    return MeshTri.init_tensor(x, y)


def fields(kind):
    def incl(x, y):
        return x ** 2 + (y - 0.5) ** 2 <= R ** 2
    if kind == "matrix":
        return (lambda x, y: MU_M + 0 * x), (lambda x, y: RHO_M + 0 * x)
    return (lambda x, y: np.where(incl(x, y), MU_I, MU_M),
            lambda x, y: np.where(incl(x, y), RHO_I, RHO_M))


def solve(W, kind, eta, n_in=161, n_out=40, ny=49):
    mesh = build(W, n_in, n_out, ny)
    vb = Basis(mesh, ElementTriP1())
    mu, rho = fields(kind)
    cm = np.sqrt(MU_M / RHO_M)
    km = eta            # period h = 1
    omega = km * cm

    @BilinearForm
    def stiff(u, v, w):
        return mu(w.x[0], w.x[1]) * dot(grad(u), grad(v))

    @BilinearForm
    def mass(u, v, w):
        return rho(w.x[0], w.x[1]) * u * v

    K = stiff.assemble(vb).astype(complex)
    M = mass.assemble(vb).astype(complex)

    fb_l = FacetBasis(mesh, ElementTriP1(), facets=mesh.facets_satisfying(lambda x: np.isclose(x[0], -W)))
    fb_r = FacetBasis(mesh, ElementTriP1(), facets=mesh.facets_satisfying(lambda x: np.isclose(x[0], W)))

    @BilinearForm
    def bmass(u, v, w):
        return u * v

    @LinearForm
    def bload(v, w):
        return v

    Ml = bmass.assemble(fb_l).astype(complex)
    Mr = bmass.assemble(fb_r).astype(complex)
    bl = bload.assemble(fb_l).astype(complex)

    A = K - omega ** 2 * M - 1j * km * MU_M * (Ml + Mr)
    b = -2j * km * MU_M * bl       # U_inc = 1 on the left facet (x=-W => e^{i k (x+W)}=1)

    P, _, _ = y_periodic_prolongation(mesh)
    Pc = P.astype(complex)
    Ar = (Pc.T @ A @ Pc).tocsc()
    br = Pc.T @ b
    Ur = spla.spsolve(Ar, br)
    U = Pc @ Ur

    p = mesh.p.T
    avg = lambda xt: U[np.abs(p[:, 0] - xt) < 1e-9].mean()
    tau = avg(W) * np.exp(-1j * km * 2 * W)
    refl = avg(-W) - 1.0
    return tau, refl, U, mesh


def main():
    print("=" * 74)
    print("  Inclusion-resolving reference: transmission across a row (normal inc.)")
    print("  mu_i/mu_m=%.2f rho_i/rho_m=%.2f R=%.4f" % (MU_I / MU_M, RHO_I / RHO_M, R))
    print("=" * 74)
    W = 6.0
    print("\n-- VERIFICATION: matrix-only must give |tau|=1, |rho|=0 --")
    for eta in (0.1, 0.3):
        tau, refl, _, _ = solve(W, "matrix", eta)
        ok = abs(abs(tau) - 1) < 1e-3 and abs(refl) < 1e-3
        print("   eta=%.2f  |tau|=%.6f  |rho|=%.2e   %s" % (eta, abs(tau), abs(refl), "PASS" if ok else "FAIL"))
    print("\n-- resolved row: transmission/reflection vs eta --")
    for eta in (0.4, 0.2, 0.1, 0.05):
        tau, refl, _, _ = solve(W, "circle", eta)
        print("   eta=%.3f  |tau|=%.6f  arg(tau)=%+.4f  |rho|=%.4e  |tau-1|=%.4e"
              % (eta, abs(tau), np.angle(tau), abs(refl), abs(tau - 1)))


if __name__ == "__main__":
    main()
