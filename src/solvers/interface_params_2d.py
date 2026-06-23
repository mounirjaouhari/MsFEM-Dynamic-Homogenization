#!/usr/bin/env python3
"""
MODEL A (single row of elastic inclusions, Marigo-faithful).  Real computation
of the effective INTERFACE jump parameters B, B2, C, C1, S from the static
boundary-layer (strip) cell problems, per MODEL_SPEC.md section 6.

Geometry: infinite strip Y_inf = R_{y1} x (0,1), 1-periodic in y2, one circular
inclusion (mu_i, rho_i) of radius R centred at (0, 1/2) in a matrix (mu_m,rho_m).
Truncated to (-L, L) x (0,1).

Correctors (MODEL_SPEC eq. in sec.6), with the SAME materials as the paper:
  normal      Q1 := u - y1, where div(mu grad u)=0, u=+-L at the x1-ends
              (matrix far field), 1-periodic in y2.
  tangential  div(mu(grad Q2 + e2))=0, Q2 1-periodic in y2, Neumann x1-ends.

Effective parameters (dimensionless, multiply the period h in the jump law):
  B  = <Q1>_{+x*} - <Q1>_{-x*}     (normal displacement-jump compliance)
  B2 = <Q2>_{+x*} - <Q2>_{-x*}     (tangential compliance)
  C  = int mu d_{y2} Q2 dy          (tangential stress/curvature coupling)
  C1 = int mu d_{y2} Q1 dy          (mixed coupling)
  S  = excess surface mass = (rho_i-rho_m)*area(incl)/period  (static)

FALSIFIABLE CHECKS (MODEL_SPEC sec.6):
  (a) centred circle  =>  B2 ~ 0 and C1 ~ 0 (reflection symmetry), C != 0, B != 0;
      energy-consistency relation B2 = -C1 (here both ~ 0).
  (b) 1D laminate (mu = mu(y1) only, full-width slab) => B2 = C = C1 = 0 exactly.
  (c) pure matrix (no inclusion) => B = B2 = C = C1 = 0.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import (MeshTri, Basis, ElementTriP1, FacetBasis,
                   BilinearForm, LinearForm, Functional)
from skfem.helpers import dot, grad

MU_M, MU_I = 12.0e9, 78.0e9          # Pa
RHO_M, RHO_I = 2500.0, 7800.0        # kg/m^3
VF = 0.25
R = np.sqrt(VF / np.pi)              # radius for 25% area fraction of the unit cell


def build_strip(L, nx, ny):
    x = np.linspace(-L, L, nx)
    y = np.linspace(0.0, 1.0, ny)
    return MeshTri.init_tensor(x, y)


def make_mu(kind="circle"):
    if kind == "circle":
        def mu(x, y):
            return np.where(x ** 2 + (y - 0.5) ** 2 <= R ** 2, MU_I, MU_M)
    elif kind == "laminate":      # mu depends on y1 only: vertical slab |x|<R, full y2 width
        def mu(x, y):
            return np.where(np.abs(x) <= R, MU_I, MU_M)
    elif kind == "tilted":        # 45-deg tilted ellipse: breaks normal/tangential symmetry
        ae, be, th = 0.28, 0.14, np.pi / 4
        c, s = np.cos(th), np.sin(th)

        def mu(x, y):
            yy = y - 0.5
            xr = c * x + s * yy
            yr = -s * x + c * yy
            return np.where((xr / ae) ** 2 + (yr / be) ** 2 <= 1.0, MU_I, MU_M)
    elif kind == "matrix":
        def mu(x, y):
            return MU_M + 0.0 * x
    return mu


def y_periodic_prolongation(mesh, tol=1e-9):
    """Identify top (y=1) nodes with bottom (y=0) nodes; return P (n x n_master)."""
    p = mesh.p.T
    n = p.shape[0]
    idx = {(round(px, 9), round(py, 9)): k for k, (px, py) in enumerate(p)}
    master = np.arange(n)
    for k in range(n):
        px, py = p[k]
        if abs(py - 1.0) < tol:
            master[k] = idx[(round(px, 9), 0.0)]
    masters = np.unique(master)
    col = {m: c for c, m in enumerate(masters)}
    rows = np.arange(n)
    cols = np.array([col[master[k]] for k in range(n)])
    P = sp.csr_matrix((np.ones(n), (rows, cols)), shape=(n, len(masters)))
    return P, masters, col


def assemble(mesh, mu):
    vb = Basis(mesh, ElementTriP1())

    @BilinearForm
    def stiff(u, v, w):
        return mu(w.x[0], w.x[1]) * dot(grad(u), grad(v))

    K = stiff.assemble(vb)
    return vb, K


def cross_mean(p, vals, xtarget, tol):
    m = np.abs(p[:, 0] - xtarget) < tol
    return vals[m].mean()


def solve_normal(mesh, vb, K, mu, L):
    """Inhomogeneous-Neumann (matrix-flux) corrector: impose the matrix normal
    flux mu_m*(e1.n) at the x-ends so the corrector u-y1 genuinely plateaus to
    its far-field constants c+- (no forced return to 0).  Then B = c+ - c-,
    measured AT the ends.  (A Dirichlet u=+-L BC under-measures B because it
    forces u-y1 -> 0 at the ends.)  C1 = int mu d_y2 (u-y1)."""
    p = mesh.p.T
    P, masters, col = y_periodic_prolongation(mesh)
    fb_r = FacetBasis(mesh, ElementTriP1(), facets=mesh.facets_satisfying(lambda x: np.isclose(x[0], L)))
    fb_l = FacetBasis(mesh, ElementTriP1(), facets=mesh.facets_satisfying(lambda x: np.isclose(x[0], -L)))
    mu_m = mu(np.array([2 * L]), np.array([0.0]))[0]    # matrix value at the ends

    @LinearForm
    def one(v, w):
        return v

    F = mu_m * (one.assemble(fb_r) - one.assemble(fb_l))   # +mu_m flux at +L, -mu_m at -L
    Kr = (P.T @ K @ P).tocsr()
    Fr = P.T @ F
    m = Kr.shape[0]
    free = np.setdiff1d(np.arange(m), [0])                 # pin one dof (mean)
    ur = np.zeros(m)
    ur[free] = spla.spsolve(Kr[free][:, free].tocsc(), Fr[free])
    u = P @ ur
    q_at = lambda xt: cross_mean(p, u - p[:, 0], xt, 1e-9)
    B = q_at(L) - q_at(-L)
    ci = vb.interpolate(u)

    @Functional
    def fy(w):
        return mu(w.x[0], w.x[1]) * grad(w["c"])[1]
    C1 = fy.assemble(vb, c=ci)
    return u, B, C1


def solve_tangential(mesh, vb, K, mu, L):
    """Corrector div(mu(grad Q2+e2))=0, y2-periodic, Neumann x-ends, pin mean."""
    P, masters, col = y_periodic_prolongation(mesh)

    @LinearForm
    def load(v, w):
        return -mu(w.x[0], w.x[1]) * grad(v)[1]     # -<mu e2.grad v>

    F = load.assemble(vb)
    Kr = (P.T @ K @ P).tocsr()
    Fr = P.T @ F
    m = Kr.shape[0]
    free = np.setdiff1d(np.arange(m), [0])          # pin one dof (mean)
    wr = np.zeros(m)
    wr[free] = spla.spsolve(Kr[free][:, free].tocsc(), Fr[free])
    W = P @ wr
    p = mesh.p.T
    xs = np.unique(np.round(p[:, 0], 9))
    xstar = xs[np.argmin(np.abs(xs - 0.6 * L))]
    q_at = lambda xt: cross_mean(p, W, xt, 1e-9)
    B2 = q_at(xstar) - q_at(-xstar)
    ci = vb.interpolate(W)

    @Functional
    def fy(w):
        return mu(w.x[0], w.x[1]) * grad(w["c"])[1]
    C = fy.assemble(vb, c=ci)
    return W, B2, C


def surface_mass():
    area_incl = np.pi * R ** 2          # per unit (y2) period
    return (RHO_I - RHO_M) * area_incl  # excess surface mass density (dimensionless cell)


def run(kind, L=4.0, nx=None, ny=129, label=""):
    nx = nx or int(round((2 * L) * (ny - 1) / 1.0)) + 1   # ~ square cells
    # cap mesh for tractability
    nx = min(nx, 8 * (ny - 1) + 1)
    mesh = build_strip(L, nx, ny)
    mu = make_mu(kind)
    vb, K = assemble(mesh, mu)
    _, B, C1 = solve_normal(mesh, vb, K, mu, L)
    _, B2, C = solve_tangential(mesh, vb, K, mu, L)
    S = surface_mass() if kind == "circle" else (surface_mass() if kind == "laminate" else 0.0)
    print(f"  [{label or kind:9s}] L={L} mesh={nx}x{ny}  "
          f"B={B:+.5f}  B2={B2:+.5e}  C={C/MU_M:+.5f}  C1={C1/MU_M:+.5e}  (C,C1 / mu_m)")
    return dict(kind=kind, L=L, nx=nx, ny=ny, B=B, B2=B2, C=C / MU_M, C1=C1 / MU_M, S=S)


if __name__ == "__main__":
    print("=" * 78)
    print("  Interface jump parameters (MODEL A, strip cell problems)  R=%.4f vf=%.3f"
          % (R, np.pi * R ** 2))
    print("=" * 78)
    print("\n-- (c) pure matrix control: all ~ 0 --")
    run("matrix", label="matrix")
    print("\n-- (b) 1D laminate control: B2 = C = C1 = 0 --")
    run("laminate", label="laminate")
    print("\n-- (a) centred circle: B,C != 0 ; B2,C1 ~ 0 (symmetry) --")
    r = run("circle", label="circle")
    print("\n-- convergence in strip half-length L (centred circle, B) --")
    for L in (2.0, 4.0, 6.0, 8.0):
        run("circle", L=L, ny=97, label=f"L={L}")
    print("\n-- (a') tilted ellipse: NONTRIVIAL test of B2 = -C1 (energy consistency) --")
    rt = run("tilted", L=4.0, ny=161, label="tilted")
    rel = abs(rt["B2"] + rt["C1"]) / max(abs(rt["B2"]), 1e-12)
    print("  B2 = %+.5f   -C1 = %+.5f   |B2+C1|/|B2| = %.3f  (-> 0 confirms B2=-C1)"
          % (rt["B2"], -rt["C1"], rel))
    print("\n  S (excess surface mass, dimensionless) = %.4f" % r["S"])
    print("  centred-circle symmetry check  B2 = -C1 :  %.2e  vs  %.2e" % (r["B2"], -r["C1"]))
