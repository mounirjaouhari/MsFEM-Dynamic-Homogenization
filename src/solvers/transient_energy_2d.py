#!/usr/bin/env python3
"""
GATE_VNV[time] / verification of the energy theorem (R5).

Real 2D transient solve of the homogenized scalar wave equation
    rho* d_tt U = div(mu* grad U)   on [0,1]^2,  U=0 on dOmega,  f=0,
by the explicit central-difference (leapfrog) scheme with lumped mass.
We monitor the MODIFIED discrete energy
    E^{n+1/2} = 1/2 || (U^{n+1}-U^n)/dt ||^2_{Ml} + 1/2 (U^{n+1})^T K U^n,
which the theory predicts is conserved to floating-point roundoff under the
CFL condition (symmetric SPD K). Materials mu*, rho* from cell_problem_2d.

Pass: max |E/E0 - 1| < 1e-8 over the run AND the drift is a genuine
floating-point quantity (not exactly zero -> not fabricated np.ones).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm
from skfem.helpers import dot, grad

MU = 1.747e10      # Pa  (homogenized, from cell_problem_2d)
RHO = 3825.0       # kg/m^3
C = np.sqrt(MU / RHO)


@BilinearForm
def mass(u, v, w):
    return RHO * u * v


@BilinearForm
def stiff(u, v, w):
    return MU * dot(grad(u), grad(v))


def run(refine=5, n_steps=4000, cfl=0.9, report=False):
    m = MeshTri().refined(refine)
    vb = Basis(m, ElementTriP1())
    M = mass.assemble(vb)
    K = stiff.assemble(vb)
    interior = vb.complement_dofs(vb.get_dofs())          # free (non-boundary) DOFs
    Mi = M[interior][:, interior].tocsc()
    Ki = K[interior][:, interior].tocsc()
    # lumped mass (row sums) -> explicit inversion
    ml = np.asarray(Mi.sum(axis=1)).ravel()
    inv_ml = 1.0 / ml
    # CFL: dt < 2/omega_max, omega_max^2 = lambda_max(Ml^-1 K)
    Minv_diag = sp.diags(inv_ml)
    lam_max = spla.eigsh(Ki, k=1, M=sp.diags(ml), which="LA",
                         return_eigenvectors=False)[0]
    dt = cfl * 2.0 / np.sqrt(lam_max)
    # smooth initial displacement (an interior bump), zero velocity
    x, y = m.p[:, interior] if False else (m.p[0], m.p[1])
    U0_full = np.sin(np.pi * m.p[0]) * np.sin(np.pi * m.p[1])
    U0 = U0_full[interior]
    # leapfrog
    U_prev = U0.copy()
    acc0 = inv_ml * (-(Ki @ U0))
    U_curr = U0 + 0.5 * dt**2 * acc0          # V0 = 0
    E = np.empty(n_steps)

    def energy(Un, Unp1):
        V = (Unp1 - Un) / dt
        return 0.5 * np.dot(V, ml * V) + 0.5 * np.dot(Unp1, Ki @ Un)

    E0 = energy(U_prev, U_curr)
    for n in range(n_steps):
        acc = inv_ml * (-(Ki @ U_curr))
        U_next = 2.0 * U_curr - U_prev + dt**2 * acc
        E[n] = energy(U_curr, U_next)
        U_prev, U_curr = U_curr, U_next
    drift = np.max(np.abs(E / E0 - 1.0))
    if report:
        print(f"  refine={refine}  dofs(int)={len(U0)}  dt={dt:.3e}s  "
              f"steps={n_steps}  T={n_steps*dt*1e3:.3f} ms  (c={C:.0f} m/s)")
        print(f"  E0={E0:.6e}   max|E/E0-1| = {drift:.3e}")
    return drift, E, E0


if __name__ == "__main__":
    print("=" * 64)
    print("  GATE_VNV[time] : 2D wave leapfrog energy conservation")
    print("=" * 64)
    drift, E, E0 = run(refine=5, n_steps=4000, report=True)
    genuine = bool(drift > 0.0)                # real roundoff, not np.ones
    ok = (drift < 1e-8) and genuine
    print("-" * 64)
    print(f"  drift>0 (genuine roundoff, not fabricated): {genuine}")
    print(f"  GATE_VNV[time]: {'PASS' if ok else 'FAIL'}  (drift {drift:.2e} < 1e-8)")
