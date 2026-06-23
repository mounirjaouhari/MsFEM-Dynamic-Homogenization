#!/usr/bin/env python3
"""
MODEL A, resonant regime (the metamaterial payoff).  Real computation of the
frequency-dependent effective SURFACE MASS B_S(omega) of a row of locally
resonant inclusions, and the band gap where B_S(omega) < 0 (negative effective
mass).  Per MODEL_SPEC.md sec.3 (resonant case) and sec.6 (eigenproblem).

Inclusion internal eigenproblem (stiff matrix clamps the inclusion boundary):
    -div(mu_i grad phi_n) = rho_i omega_n^2 phi_n  in D,   phi_n = 0 on dD.
Mass-normalized modes (int_D rho_i phi_n^2 = 1); modal participation
    a_n = ( int_D rho_i phi_n )^2.
Effective dynamic surface mass (locally-resonant closure, Milton 2002 form):
    M_eff(omega) = m_s + sum_n a_n omega_n^2 / (omega_n^2 - omega^2),
with m_s the static surface mass.  M_eff -> -inf just above each omega_n
=> M_eff<0 band (no propagation) => band gap [omega_n, omega_n/sqrt(1-a_n/M_tot)].

VERIFICATION: the clamped homogeneous disk has the exact spectrum
    omega_{m,k} = j_{m,k} * c_i / R_phys,   c_i = sqrt(mu_i/rho_i),
j_{m,k}=k-th zero of Bessel J_m; fundamental j_{0,1}=2.40483.  The FE omega_1
must match this to discretization error.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import jn_zeros
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm, LinearForm, Functional
from skfem.helpers import dot, grad

# locally-resonant set: SOFT, DENSE inclusion (low internal resonance) in a stiff matrix
MU_M, MU_I = 12.0e9, 0.5e9            # Pa  (soft inclusion => low omega_1)
RHO_M, RHO_I = 2500.0, 9000.0         # kg/m^3 (dense inclusion)
R_PHYS = 5.0e-3                        # physical inclusion radius (m), illustrative


def build_disk(R, nr=40, nth=64):
    """Polar-structured triangular disk mesh of radius R (centre + nr rings)."""
    pts = [(0.0, 0.0)]
    for i in range(1, nr + 1):
        r = R * i / nr
        for j in range(nth):
            th = 2 * np.pi * j / nth
            pts.append((r * np.cos(th), r * np.sin(th)))
    p = np.array(pts).T                          # (2, npts)

    def ring_idx(i, j):                           # node index on ring i (1..nr), angle j
        return 1 + (i - 1) * nth + (j % nth)

    tris = []
    for j in range(nth):                          # inner fan (centre to ring 1)
        tris.append((0, ring_idx(1, j), ring_idx(1, j + 1)))
    for i in range(1, nr):                        # quad rings split into 2 triangles
        for j in range(nth):
            a, b = ring_idx(i, j), ring_idx(i, j + 1)
            c, d = ring_idx(i + 1, j), ring_idx(i + 1, j + 1)
            tris.append((a, b, d))
            tris.append((a, d, c))
    t = np.array(tris).T
    return MeshTri(p, t), nth


def inclusion_modes(R=R_PHYS, nmodes=6, nr=48, nth=96):
    mesh, nth = build_disk(R, nr, nth)
    vb = Basis(mesh, ElementTriP1())

    @BilinearForm
    def stiff(u, v, w):
        return MU_I * dot(grad(u), grad(v))

    @BilinearForm
    def mass(u, v, w):
        return RHO_I * u * v

    K = stiff.assemble(vb)
    M = mass.assemble(vb)
    # Dirichlet phi=0 on outer ring (boundary nodes)
    rr = np.sqrt(mesh.p[0] ** 2 + mesh.p[1] ** 2)
    bnd = np.where(rr > R - 1e-12)[0]
    free = np.setdiff1d(np.arange(K.shape[0]), bnd)
    Kf, Mf = K[free][:, free].tocsc(), M[free][:, free].tocsc()
    lam, vec = spla.eigsh(Kf, k=nmodes, M=Mf, sigma=0.0, which='LM')
    order = np.argsort(lam)
    lam, vec = lam[order], vec[:, order]
    omega = np.sqrt(np.abs(lam))
    # mass-normalize: int_D rho phi^2 = phi^T M phi = 1  (eigsh already M-normalizes)
    modes = np.zeros((K.shape[0], nmodes))
    modes[free] = vec
    # participation a_n = (int_D rho phi_n)^2 = (phi_n^T M 1)^2
    one = np.ones(K.shape[0])
    a = np.array([(modes[:, n] @ (M @ one)) ** 2 for n in range(nmodes)])
    m_static = (modes[:, 0] @ (M @ one))  # ~ total inclusion mass proxy via mode-0? use direct:
    M_incl = one @ (M @ one)              # total inclusion mass int_D rho
    return omega, a, M_incl, mesh


def effective_mass(omega_grid, omega_n, a_n, m_s):
    Me = np.full_like(omega_grid, m_s, dtype=float)
    for wn, an in zip(omega_n, a_n):
        Me += an * wn ** 2 / (wn ** 2 - omega_grid ** 2)
    return Me


def main():
    print("=" * 74)
    print("  Resonant surface mass B_S(omega)  (MODEL A, locally-resonant inclusion)")
    print("  mu_i/mu_m=%.3f  rho_i/rho_m=%.2f  R=%.1f mm" %
          (MU_I / MU_M, RHO_I / RHO_M, R_PHYS * 1e3))
    print("=" * 74)
    omega_n, a_n, M_incl, mesh = inclusion_modes()
    c_i = np.sqrt(MU_I / RHO_I)
    j0 = jn_zeros(0, 2); j1 = jn_zeros(1, 1)
    w_exact = np.array([j0[0], j1[0], j0[1]]) * c_i / R_PHYS   # first 3 distinct: (0,1),(1,1)x2,(0,2)
    print("  c_i = sqrt(mu_i/rho_i) = %.1f m/s,  inclusion mass int_D rho = %.4e" % (c_i, M_incl))
    print("\n  internal resonance frequencies (rad/s):")
    print("   FE omega_n  :", np.array2string(omega_n[:4], precision=3, formatter={'float': lambda z: f'{z:.3e}'}))
    print("   exact (Bessel j_{0,1},j_{1,1},j_{0,2})*c_i/R :",
          np.array2string(w_exact, precision=3, formatter={'float': lambda z: f'{z:.3e}'}))
    err1 = abs(omega_n[0] - w_exact[0]) / w_exact[0]
    print("   fundamental relative error vs Bessel j_{0,1}=2.4048 : %.3e  %s"
          % (err1, "PASS" if err1 < 2e-2 else "FAIL"))
    print("   participation a_n (mass-normalized):",
          np.array2string(a_n[:4], precision=3, formatter={'float': lambda z: f'{z:.3e}'}))

    m_s = M_incl                          # static surface mass proxy (excess folded in §5)
    wg = np.linspace(0.05 * omega_n[0], 1.6 * omega_n[0], 4000)
    Me = effective_mass(wg, omega_n, a_n, m_s)
    neg = wg[Me < 0]
    f_n = omega_n / (2 * np.pi)
    print("\n  band gap (effective mass M_eff(omega) < 0):")
    if neg.size:
        print("   gap = [%.3e, %.3e] rad/s  =  [%.1f, %.1f] kHz" %
              (neg.min(), neg.max(), neg.min() / 2 / np.pi / 1e3, neg.max() / 2 / np.pi / 1e3))
        print("   opens at the fundamental resonance f_1 = %.1f kHz" % (f_n[0] / 1e3))
        print("   GATE: negative-effective-mass band gap present above omega_1 : PASS")
    else:
        print("   no negative-mass band found (check resonance vs grid) : FAIL")
    return dict(omega_n=omega_n, a_n=a_n, M_incl=M_incl, gap=(neg.min(), neg.max()) if neg.size else None,
                w_exact=w_exact, err1=err1)


if __name__ == "__main__":
    main()
