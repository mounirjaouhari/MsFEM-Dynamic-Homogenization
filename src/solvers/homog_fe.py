#!/usr/bin/env python3
"""
REAL finite-element solve of the HOMOGENIZED model (Model A): matrix bulk on
both sides of Gamma={x=0}, coupled by the symmetric Ventcel/spring interface that
imposes the Marigo displacement jump  [U] = h*B <d_x U>.  No inclusion is meshed;
the row is represented solely by the interface condition.  This is the paper's
actual MsFEM--Ventcel method, solved (not drawn).

Static unit-gradient (corrector) version: two matrix half-strips (-L,0)x(0,1) and
(0,L)x(0,1), y-periodic, driven by the matrix normal flux +-mu_m at x=+-L; the
interface couples them by  +kappa * int_Gamma [U][V],  kappa = +mu_m/(h*B)
(<0 for B<0, the Option-A negative surface stiffness).  Verification: the
produced jump <U_R - U_L>_Gamma must equal h*B.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import (MeshTri, Basis, ElementTriP1, FacetBasis, BilinearForm, LinearForm)
from skfem.helpers import dot, grad
from interface_params_2d import y_periodic_prolongation, MU_M

H = 1.0   # period (nondim)


def half_strip(x0, x1, nx, ny):
    return MeshTri.init_tensor(np.linspace(x0, x1, nx), np.linspace(0, 1, ny))


def assemble_side(mesh, flux_x, sign):
    """Matrix stiffness + matrix-flux Neumann load at the outer end x=flux_x."""
    vb = Basis(mesh, ElementTriP1())

    @BilinearForm
    def stiff(u, v, w):
        return MU_M * dot(grad(u), grad(v))

    K = stiff.assemble(vb)
    fb = FacetBasis(mesh, ElementTriP1(), facets=mesh.facets_satisfying(lambda x: np.isclose(x[0], flux_x)))

    @LinearForm
    def one(v, w):
        return v

    F = sign * MU_M * one.assemble(fb)        # +mu_m flux at right end, -mu_m at left end
    return vb, K, F


def iface_dofs(mesh, P, col, xface):
    """Reduced-DOF indices of the x=xface nodes, sorted by y, with their y and dy (lumped)."""
    p = mesh.p.T
    nodes = np.where(np.abs(p[:, 0] - xface) < 1e-9)[0]
    red = P[nodes].nonzero()[1]               # reduced col per node
    ys = p[nodes, 1]
    order = np.argsort(ys)
    nodes, red, ys = nodes[order], red[order], ys[order]
    # lumped 1D mass (nodal length), periodic in y
    dy = np.gradient(ys)
    return red, ys, dy


def solve_homog_static(B, L=4.0, nx=120, ny=61):
    kappa = MU_M / (H * B)                     # spring s.t. [U]=Sigma1/kappa = h*B under unit flux
    mL = half_strip(-L, 0.0, nx, ny)
    mR = half_strip(0.0, L, nx, ny)
    vbL, KL, FL = assemble_side(mL, -L, -1.0)
    vbR, KR, FR = assemble_side(mR, +L, +1.0)
    PL, mastersL, colL = y_periodic_prolongation(mL)
    PR, mastersR, colR = y_periodic_prolongation(mR)
    KLr = (PL.T @ KL @ PL).tocsr(); FLr = PL.T @ FL
    KRr = (PR.T @ KR @ PR).tocsr(); FRr = PR.T @ FR
    nL, nR = KLr.shape[0], KRr.shape[0]
    # interface dofs (left mesh x=0 ; right mesh x=0), matched by y
    redL, yL, dyL = iface_dofs(mL, PL, colL, 0.0)
    redR, yR, dyR = iface_dofs(mR, PR, colR, 0.0)
    assert np.allclose(yL, yR)
    dy = 0.5 * (dyL + dyR)
    # global block system  [KLr , 0 ; 0 , KRr] + interface coupling kappa*[U][V]
    A = sp.bmat([[KLr, None], [None, KRr]], format="lil")
    F = np.concatenate([FLr, FRr])
    for a, b, d in zip(redL, redR + nL, dy):   # [U]=U_R-U_L ; coupling kappa*d*([..])
        A[a, a] += kappa * d; A[b, b] += kappa * d
        A[a, b] += -kappa * d; A[b, a] += -kappa * d
    A = A.tocsr()
    free = np.setdiff1d(np.arange(nL + nR), [0])   # pin one dof (pure Neumann)
    u = np.zeros(nL + nR)
    u[free] = spla.spsolve(A[free][:, free].tocsc(), F[free])
    UL = PL @ u[:nL]; UR = PR @ u[nL:]
    # interface nodes in FULL node indexing (x=0), sorted by y, matched
    pL = mL.p.T; pR = mR.p.T
    iL = np.where(np.abs(pL[:, 0]) < 1e-9)[0]; iL = iL[np.argsort(pL[iL, 1])]
    iR = np.where(np.abs(pR[:, 0]) < 1e-9)[0]; iR = iR[np.argsort(pR[iR, 1])]
    jump = float((UR[iR] - UL[iL]).mean())
    return dict(mL=mL, mR=mR, UL=UL, UR=UR, jump=jump, B=B)


def main():
    from interface_params_2d import build_strip, make_mu, assemble, solve_normal
    # get the cell B (circle, contrast 6.5)
    mesh = build_strip(4.0, 481, 61); mu = make_mu("circle")
    vb, K = assemble(mesh, mu); _, B, _ = solve_normal(mesh, vb, K, mu, 4.0)
    print("cell B (circle) = %.4f" % B)
    r = solve_homog_static(B)
    print("homogenized FE solve: produced jump <U_R-U_L> = %.4f   (target h*B = %.4f)"
          % (r["jump"], H * B))
    err = abs(r["jump"] - H * B) / abs(H * B)
    print("relative error = %.2e   %s" % (err, "PASS" if err < 0.05 else "CHECK"))


if __name__ == "__main__":
    main()
