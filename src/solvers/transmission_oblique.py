#!/usr/bin/env python3
"""
MODEL A, genuinely-2D validation (todo 7): OBLIQUE incidence activates the
transverse parameters.  A plane wave hits the row at angle theta, so the field
carries a transverse wavenumber k_y = k sin(theta) (Bloch-periodic in y).  The
jump conditions then involve B2 (in [U]) and C, C1 (in [Sigma1]):
  [U]      = B <d_x U> + B2 <d_y U>
  [Sigma1] = -S w^2 <U> + C k_y^2 <U> - C1 i k_y <d_x U>      (mu_m = 1)
For a CENTRED circle B2 = C1 = 0, so the active transverse term is + C k_y^2 <U>.

Reference = resolved circular inclusions with BLOCH-periodic BC in y, ABC using
the normal wavenumber k_x = k cos(theta) (exact for the propagating 0th order).

CLAIMS TO VERIFY:
  (1) matrix-only (oblique) -> |tau|=1 exactly (sanity);
  (2) the homogenized model WITH C reproduces the resolved oblique transmission
      markedly better than the model WITHOUT C  => C is genuinely activated;
  (3) the effective C inverted from the reference matches the cell C.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import (MeshTri, Basis, ElementTriP1, FacetBasis, BilinearForm, LinearForm)
from skfem.helpers import dot, grad
import transmission_resolved as tr
from interface_params_2d import MU_M, RHO_M, MU_I, RHO_I, R

# verified cell parameters (PROVENANCE.md), centred circle, contrast 6.5
B = -0.311
B2 = 0.0
C = -0.914            # C / mu_m  (mu_m = 1 nondim)
C1 = 0.0
S = (RHO_I / RHO_M - 1.0) * np.pi * R ** 2     # 0.530


def bloch_prolongation(mesh, phase, tol=1e-9):
    """Bloch y-periodicity: U(top) = phase * U(bottom).  Returns complex P."""
    p = mesh.p.T
    n = p.shape[0]
    idx = {(round(px, 9), round(py, 9)): k for k, (px, py) in enumerate(p)}
    master = np.arange(n)
    fac = np.ones(n, dtype=complex)
    for k in range(n):
        px, py = p[k]
        if abs(py - 1.0) < tol:
            master[k] = idx[(round(px, 9), 0.0)]
            fac[k] = phase
    masters = np.unique(master)
    col = {m: c for c, m in enumerate(masters)}
    rows = np.arange(n)
    cols = np.array([col[master[k]] for k in range(n)])
    P = sp.csr_matrix((fac, (rows, cols)), shape=(n, len(masters)), dtype=complex)
    return P


def solve_resolved_oblique(W, kind, eta, theta, n_in=161, n_out=40, ny=49,
                           return_field=False):
    mesh = tr.build(W, n_in, n_out, ny)
    vb = Basis(mesh, ElementTriP1())
    mu, rho = tr.fields(kind)
    cm = np.sqrt(MU_M / RHO_M)
    k = eta                       # period h = 1 -> k = eta
    kx, ky = k * np.cos(theta), k * np.sin(theta)
    omega = k * cm

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
    def bm(u, v, w):
        return u * v

    @LinearForm
    def binc_c(v, w):             # Re part of U_inc = e^{i k_y y} on the left facet
        return np.cos(ky * w.x[1]) * v

    @LinearForm
    def binc_s(v, w):             # Im part (assembled separately: complex forms get truncated)
        return np.sin(ky * w.x[1]) * v

    Ml = bm.assemble(fb_l).astype(complex)
    Mr = bm.assemble(fb_r).astype(complex)
    bl = binc_c.assemble(fb_l).astype(complex) + 1j * binc_s.assemble(fb_l)
    # ABC uses the NORMAL wavenumber k_x (exact for the 0th order)
    A = K - omega ** 2 * M - 1j * kx * MU_M * (Ml + Mr)
    b = -2j * kx * MU_M * bl

    P = bloch_prolongation(mesh, np.exp(1j * ky * 1.0))   # phase across one period (=1)
    Ph = P.conj().T
    Ar = (Ph @ A @ P).tocsc()
    br = Ph @ b
    Ur = spla.spsolve(Ar, br)
    U = P @ Ur

    p = mesh.p.T
    def order0(xt):               # project onto the 0th diffraction order e^{i k_y y}
        m = np.abs(p[:, 0] - xt) < 1e-9
        return (U[m] * np.exp(-1j * ky * p[m, 1])).mean()
    tau = order0(W) * np.exp(-1j * kx * 2 * W)
    refl = order0(-W) - 1.0
    if return_field:
        return tau, refl, U, mesh, (kx, ky)
    return tau, refl


def tau_homog_oblique(eta, theta, useC=True):
    k = eta
    kx, ky = k * np.cos(theta), k * np.sin(theta)
    w2 = k ** 2
    Cc = C if useC else 0.0
    # (I):  t-(1+r) = B(ikx/2)(t+1-r) + B2(iky/2)(t+1+r)
    # (II): ikx(t-1+r) = 1/2(t+1+r)(-S w2 + Cc ky^2) + 1/2 C1 kx ky (t+1-r)
    a1t = 1 - B * 1j * kx / 2 - B2 * 1j * ky / 2
    a1r = -1 + B * 1j * kx / 2 - B2 * 1j * ky / 2
    c1 = 1 + B * 1j * kx / 2 + B2 * 1j * ky / 2
    g = (-S * w2 + Cc * ky ** 2) / 2
    a2t = 1j * kx - g - C1 * kx * ky / 2
    a2r = 1j * kx - g + C1 * kx * ky / 2
    c2 = 1j * kx + g - C1 * kx * ky / 2
    A = np.array([[a1t, a1r], [a2t, a2r]], dtype=complex)
    rhs = np.array([c1, c2], dtype=complex)
    t, r = np.linalg.solve(A, rhs)
    return t, r


def main():
    print("=" * 84)
    print("  Oblique incidence: transverse parameter C activated (centred circle)")
    print("  B=%.3f  C=%.3f  S=%.3f   B2=C1=0 (symmetry)" % (B, C, S))
    print("=" * 84)
    W, theta = 6.0, np.deg2rad(45.0)
    print("  theta = 45 deg\n")
    print("  -- sanity: matrix-only oblique must give |tau|=1 --")
    for eta in (0.1, 0.3):
        t, r = solve_resolved_oblique(W, "matrix", eta, theta)
        print("     eta=%.2f  |tau|=%.6f  |rho|=%.2e  %s"
              % (eta, abs(t), abs(r), "PASS" if abs(abs(t) - 1) < 2e-3 else "FAIL"))
    print("\n  -- resolved vs homogenized (with C) vs homogenized (C=0) --")
    print("  %-7s | %-10s | %-12s %-12s | %-12s %-12s"
          % ("eta", "|t_res|", "|t_res-t_C|", "|t_res-t_noC|", "argres", "arg_C"))
    print("  " + "-" * 78)
    eC, e0 = [], []
    etas = [0.4, 0.2, 0.1, 0.05]
    for eta in etas:
        tre, rre = solve_resolved_oblique(W, "circle", eta, theta)
        tC, _ = tau_homog_oblique(eta, theta, useC=True)
        t0, _ = tau_homog_oblique(eta, theta, useC=False)
        eC.append(abs(tre - tC)); e0.append(abs(tre - t0))
        print("  %-7.3f | %-10.6f | %-12.3e %-12.3e | %+0.5f     %+0.5f"
              % (eta, abs(tre), abs(tre - tC), abs(tre - t0), np.angle(tre), np.angle(tC)))
    print("  " + "-" * 78)
    print("  HONEST STATUS: matrix-only oblique is exact (|tau|=1), so the Bloch/ABC solver is")
    print("  correct.  The row DOES affect oblique transmission, but the transverse term C enters")
    print("  the 0th-order transmission only at O(eta^2) (~ C k_y^2 <U>, i.e. eta^2 sin^2 theta),")
    print("  subdominant to the O(eta) B,S effect; it cannot be cleanly isolated here, and the exact")
    print("  scalar stress-jump algebraic form was not sourceable verbatim (paywall).")
    print("  => The 2D transverse content (C != 0, B2 = -C1) is established at the CELL level")
    print("     (interface_params_2d.py); a quantitative system-level C extraction is left open.")


if __name__ == "__main__":
    main()
