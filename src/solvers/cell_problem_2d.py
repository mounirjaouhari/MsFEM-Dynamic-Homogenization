#!/usr/bin/env python3
"""
2D periodic cell problem on Y=[0,1]^2 with a circular inclusion (antiplane
shear / 2D conductivity).  Solves the corrector problems

    div_y( mu(y) ( grad chi^j + e_j ) ) = 0 ,  chi^j Y-periodic, <chi^j>=0,

and forms the homogenized tensor  mu*_{ij} = < mu ( d_i chi^j + delta_ij ) >.

Verification (GATE_VNV, cell part): mu* must be SPD, ~isotropic for a centred
circle (mu*_11 ~ mu*_22, mu*_12 ~ 0), and lie between the Voigt and Reuss means
and inside the 2D Hashin-Shtrikman bounds. Materials = paper's test case
(epoxy matrix mu_m=12 GPa, aluminium inclusion mu_i=78 GPa, vf=25%).

Real computation; periodicity imposed by an explicit DOF reduction.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm, LinearForm, Functional
from skfem.helpers import dot, grad

MU_M, MU_I = 12.0e9, 78.0e9          # Pa
VF = 0.25
R = np.sqrt(VF / np.pi)              # inclusion radius for the target vol. frac.
CX = CY = 0.5


def mu_at(x, y):
    return np.where((x - CX) ** 2 + (y - CY) ** 2 <= R ** 2, MU_I, MU_M)


@BilinearForm
def stiff_mu(u, v, w):
    x, y = w.x
    return mu_at(x, y) * dot(grad(u), grad(v))


def load_dir(j):
    @LinearForm
    def f(v, w):
        x, y = w.x
        g = grad(v)
        return -mu_at(x, y) * g[j]      # -< mu e_j . grad v >
    return f


def mu_times_dchi(i):
    @Functional
    def F(w):
        x, y = w.x
        g = grad(w["chi"])
        return mu_at(x, y) * g[i]
    return F


@Functional
def mu_avg(w):
    x, y = w.x
    return mu_at(x, y) + 0.0 * w["chi"]   # keep shape via dummy field use


def periodic_prolongation(mesh, tol=1e-9):
    """Map every node to its periodic master (right->left, top->bottom,
    corners->origin) and return sparse P (n_full x n_master)."""
    p = mesh.p.T                                   # (n,2)
    n = p.shape[0]
    coord_to_idx = {}
    for k in range(n):
        coord_to_idx[(round(p[k, 0], 9), round(p[k, 1], 9))] = k
    master = np.empty(n, dtype=int)
    for k in range(n):
        mx, my = p[k, 0], p[k, 1]
        if abs(mx - 1.0) < tol:
            mx = 0.0
        if abs(my - 1.0) < tol:
            my = 0.0
        master[k] = coord_to_idx[(round(mx, 9), round(my, 9))]
    masters = np.unique(master)
    col = {m: c for c, m in enumerate(masters)}
    rows = np.arange(n)
    cols = np.array([col[master[k]] for k in range(n)])
    P = sp.csr_matrix((np.ones(n), (rows, cols)), shape=(n, len(masters)))
    return P, masters, col


def solve_corrector(basis, K, P, fvec, pin_red=0):
    """Reduce to periodic DOFs, pin one (mean fixed later is irrelevant for mu*)."""
    Kr = (P.T @ K @ P).tocsr()
    Fr = P.T @ fvec
    m = Kr.shape[0]
    free = np.setdiff1d(np.arange(m), [pin_red])
    chi_red = np.zeros(m)
    chi_red[free] = spla.spsolve(Kr[free][:, free].tocsc(), Fr[free])
    return P @ chi_red                              # full nodal vector


def main():
    print("=" * 64)
    print("  2D cell problem  (R=%.4f, vf=%.3f, mu_i/mu_m=%.2f)"
          % (R, np.pi * R ** 2, MU_I / MU_M))
    print("=" * 64)
    mesh = MeshTri().refined(7)                     # h=1/128
    vb = Basis(mesh, ElementTriP1())
    K = stiff_mu.assemble(vb)
    P, masters, col = periodic_prolongation(mesh)

    chi = [solve_corrector(vb, K, P, load_dir(j).assemble(vb)) for j in (0, 1)]

    mu_bar = mu_avg.assemble(vb, chi=vb.interpolate(chi[0]))   # <mu>
    mustar = np.zeros((2, 2))
    for j in (0, 1):
        ci = vb.interpolate(chi[j])
        for i in (0, 1):
            mustar[i, j] = mu_times_dchi(i).assemble(vb, chi=ci) + (mu_bar if i == j else 0.0)

    # closed-form bounds
    voigt = (1 - VF) * MU_M + VF * MU_I
    reuss = 1.0 / ((1 - VF) / MU_M + VF / MU_I)
    # 2D Hashin-Shtrikman (d=2): denominator uses 2*mu of the reference phase
    hs_lo = MU_M + VF / (1.0 / (MU_I - MU_M) + (1 - VF) / (2 * MU_M))
    hs_hi = MU_I + (1 - VF) / (1.0 / (MU_M - MU_I) + VF / (2 * MU_I))
    eig = np.linalg.eigvalsh(0.5 * (mustar + mustar.T))

    print("  <mu>           = %.4e Pa" % mu_bar)
    print("  mu* =\n   [[% .4e, % .4e],\n    [% .4e, % .4e]] Pa"
          % (mustar[0, 0], mustar[0, 1], mustar[1, 0], mustar[1, 1]))
    print("  eigenvalues    = %.4e , %.4e Pa" % (eig[0], eig[1]))
    print("  Reuss  = %.4e   HS- = %.4e   HS+ = %.4e   Voigt = %.4e"
          % (reuss, hs_lo, hs_hi, voigt))
    spd = bool(eig[0] > 0)
    iso = abs(mustar[0, 0] - mustar[1, 1]) / mustar[0, 0] < 0.02 and abs(mustar[0, 1]) / mustar[0, 0] < 0.02
    in_vr = bool(reuss - 1 <= eig[0] and eig[1] <= voigt + 1)
    in_hs = bool(hs_lo * 0.999 <= eig[0] and eig[1] <= hs_hi * 1.001)
    print("-" * 64)
    print("  SPD: %s | isotropic(centred circle): %s | within Voigt-Reuss: %s | within HS: %s"
          % (spd, iso, in_vr, in_hs))
    print("  GATE_VNV[cell]: %s" % ("PASS" if (spd and iso and in_vr and in_hs) else "FAIL"))
    return mustar, mu_bar, eig, (reuss, hs_lo, hs_hi, voigt)


if __name__ == "__main__":
    main()
