#!/usr/bin/env python3
"""
Transient scattering validation — Figure 6b (fig_scattering_3times.pdf).

DNS: full 2D FEM simulation on a tensor mesh with the circular inclusion
     resolved through element-wise material assignment.
     Leapfrog time integration (row-lumped mass), y-periodic Neumann BCs.

Homogenized (EI): TWO separate half-domains (x < -R and x > R) coupled at
     x = ±R by discrete VENTCEL jump conditions (Enlarged-Interface model,
     e=2R).  Parameters B=-0.31, S=0.53 as in the article (Table 2, Model A).

Visual layout (3 x 3):
  Row 0 — DNS 2D field + circle outlines       at t1, t2, t3
  Row 1 — EI homogenized field + slab lines    at t1, t2, t3
  Row 2 — 1D cross-section at y=H/2 (DNS vs EI)  at t1, t2, t3

Parameters (structural regime):
  mu_i / mu_m = 6.5   (MU_I=78e9, MU_M=12e9 Pa)
  rho_i / rho_m = 3.1  (RHO_I=7800, RHO_M=2500 kg/m³)
  R = sqrt(VF/pi) ≈ 0.282  (VF=0.25)
  B = -0.31,  S = 0.53  (interface jump parameters, dimensionless)
  e = 2R  (Enlarged-Interface thickness)

Output: figures/fig_scattering_3times.pdf / .png
"""
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _d in ('solvers', 'utils'):
    _p = _os.path.join(_SRC, _d)
    if _p not in _sys.path: _sys.path.insert(0, _p)
del _sys, _os, _SRC, _d, _p

import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import pubstyle; pubstyle.apply()
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm
from skfem.helpers import dot, grad
from interface_params_2d import MU_M, MU_I, RHO_M, RHO_I, R

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)

# ── physical parameters ───────────────────────────────────────────────────────
L        = 3.0       # half-domain in x
H_period = 1.0
C_M      = np.sqrt(MU_M / RHO_M)

B   = -0.31    # normal displacement-jump compliance (dimensionless, Table 2)
S   =  0.53    # excess surface mass (dimensionless, Table 2)

def incl(x, y):
    return x**2 + (y - 0.5)**2 <= R**2

# ── DNS: 2D FEM tensor mesh, inclusion resolved ───────────────────────────────
def get_dns_matrices():
    nx, ny = 201, 31
    m = MeshTri.init_tensor(np.linspace(-L, L, nx), np.linspace(0.0, 1.0, ny))
    vb = Basis(m, ElementTriP1())

    @BilinearForm
    def stiff(u, v, w):
        return np.where(incl(w.x[0], w.x[1]), MU_I, MU_M) * dot(grad(u), grad(v))

    @BilinearForm
    def mass(u, v, w):
        return np.where(incl(w.x[0], w.x[1]), RHO_I, RHO_M) * u * v

    M = mass.assemble(vb)
    K = stiff.assemble(vb)
    return m, K, M


# ── Homogenized: two half-domains + Ventcel jump conditions at x=±R ───────────
def get_homog_matrices():
    e_gap = 2.0 * R
    nx, ny = 120, 31
    mL = MeshTri.init_tensor(np.linspace(-L, -e_gap / 2, nx // 2), np.linspace(0.0, 1.0, ny))
    mR = MeshTri.init_tensor(np.linspace( e_gap / 2,  L, nx // 2), np.linspace(0.0, 1.0, ny))

    vbL = Basis(mL, ElementTriP1())
    vbR = Basis(mR, ElementTriP1())

    @BilinearForm
    def stiff_hom(u, v, w):
        return MU_M * dot(grad(u), grad(v))

    @BilinearForm
    def mass_hom(u, v, w):
        return RHO_M * u * v

    KL = stiff_hom.assemble(vbL); ML = mass_hom.assemble(vbL)
    KR = stiff_hom.assemble(vbR); MR = mass_hom.assemble(vbR)
    nL = KL.shape[0]

    def ifaceL(mesh):
        p = mesh.p.T
        n = np.where(np.abs(p[:, 0] - (-e_gap / 2)) < 1e-9)[0]
        ys = p[n, 1]; o = np.argsort(ys)
        return n[o], ys[o], np.gradient(ys[o])

    def ifaceR(mesh):
        p = mesh.p.T
        n = np.where(np.abs(p[:, 0] - (e_gap / 2)) < 1e-9)[0]
        ys = p[n, 1]; o = np.argsort(ys)
        return n[o], ys[o], np.gradient(ys[o])

    rL, yL, dL = ifaceL(mL)
    rR, yR, dR = ifaceR(mR)
    dy = 0.5 * (dL + dR)

    K_hom = sp.bmat([[KL, None], [None, KR]], format="lil")
    M_hom = sp.bmat([[ML, None], [None, MR]], format="lil")

    B_enl = B + e_gap / H_period
    S_enl = S + e_gap / H_period
    kap   = MU_M / (H_period * B_enl)
    sm    = H_period * S_enl * RHO_M

    for a, b, d in zip(rL, rR + nL, dy):
        K_hom[a, a] += kap * d;  K_hom[b, b] += kap * d
        K_hom[a, b] -= kap * d;  K_hom[b, a] -= kap * d
        M_hom[a, a] += 0.25 * sm * d;  M_hom[b, b] += 0.25 * sm * d
        M_hom[a, b] += 0.25 * sm * d;  M_hom[b, a] += 0.25 * sm * d

    class HomogMesh:
        def __init__(self, mL, mR):
            self.mL = mL; self.mR = mR

    return HomogMesh(mL, mR), K_hom.tocsc(), M_hom.tocsc(), nL


# ── leapfrog solver ───────────────────────────────────────────────────────────
def _leapfrog(K, M, x_coords, save_times):
    """Row-lumped explicit leapfrog.  Returns list of U snapshots at save_times."""
    ml    = np.asarray(M.sum(axis=1)).ravel()
    inv_ml = 1.0 / ml
    lam   = spla.eigsh(K.tocsc(), k=1, M=sp.diags(ml), which="LA",
                       return_eigenvectors=False)[0]
    dt    = 0.9 * 2.0 / np.sqrt(lam)

    x0, sigma = -1.5, 0.25
    U0 = np.exp(-(x_coords - x0)**2 / (2.0 * sigma**2))
    V0 = C_M * (x_coords - x0) / (sigma**2) * U0

    acc0   = inv_ml * (-(K @ U0))
    U_curr = U0.copy()
    U_prev = U0 - dt * V0 + 0.5 * dt**2 * acc0

    frames = []
    t = 0.0
    for target_t in save_times:
        while t < target_t:
            acc    = inv_ml * (-(K @ U_curr))
            U_next = 2.0 * U_curr - U_prev + dt**2 * acc
            U_prev, U_curr = U_curr, U_next
            t += dt
        frames.append(U_curr.copy())
    return frames


# ── triangulations for multi-period display ───────────────────────────────────
def _tiled_tri(mesh, nrep=3):
    p = mesh.p; base_t = mesh.t.T
    xs, ys, tris = [], [], []
    for n in range(nrep):
        xs.append(p[0]); ys.append(p[1] + n)
        tris.append(base_t + n * p.shape[1])
    return mtri.Triangulation(
        np.concatenate(xs), np.concatenate(ys),
        np.concatenate(tris, axis=0)
    )


# ── figure ────────────────────────────────────────────────────────────────────
def make_figure():
    save_times = [450e-6, 850e-6, 1400e-6]

    print("  Assembling DNS matrices...")
    m_dns, K_dns, M_dns = get_dns_matrices()

    print("  Assembling homogenized matrices (Ventcel jump conditions)...")
    m_hom, K_hom, M_hom, nL = get_homog_matrices()

    print("  Running DNS leapfrog...")
    frames_dns = _leapfrog(K_dns, M_dns, m_dns.p[0], save_times)

    print("  Running EI homogenized leapfrog...")
    p_hom_all = np.concatenate([m_hom.mL.p[0], m_hom.mR.p[0]])
    frames_hom = _leapfrog(K_hom, M_hom, p_hom_all, save_times)

    # ── triangulations for 3-period tiling ───────────────────────────────────
    nrep    = 3
    tri_dns = _tiled_tri(m_dns, nrep)
    tri_hL  = _tiled_tri(m_hom.mL, nrep)
    tri_hR  = _tiled_tri(m_hom.mR, nrep)

    # ── 1D slice nodes at y=0.5 ───────────────────────────────────────────────
    p_dns = m_dns.p.T
    idx_dns = np.where(np.abs(p_dns[:, 1] - 0.5) < 1e-6)[0]
    idx_dns = idx_dns[np.argsort(p_dns[idx_dns, 0])]
    x_dns   = p_dns[idx_dns, 0]

    pL = m_hom.mL.p.T; pR = m_hom.mR.p.T
    idxL = np.where(np.abs(pL[:, 1] - 0.5) < 1e-6)[0]
    idxL = idxL[np.argsort(pL[idxL, 0])]; xL = pL[idxL, 0]
    idxR = np.where(np.abs(pR[:, 1] - 0.5) < 1e-6)[0]
    idxR = idxR[np.argsort(pR[idxR, 0])]; xR = pR[idxR, 0]

    # ── plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(
        3, 3, figsize=(15, 10),
        gridspec_kw={'height_ratios': [1, 1, 1.2]}
    )
    fig.subplots_adjust(hspace=0.30, wspace=0.15)

    lim = 0.4
    lvl = np.linspace(-lim, lim, 41)
    th  = np.linspace(0.0, 2.0 * np.pi, 60)

    for j, (target_t, U_dns, U_hom) in enumerate(
            zip(save_times, frames_dns, frames_hom)):
        ax0 = axes[0, j]; ax1 = axes[1, j]; ax2 = axes[2, j]

        # ── Row 0: DNS ────────────────────────────────────────────────────────
        V_dns = np.concatenate([U_dns for _ in range(nrep)])
        ax0.tricontourf(tri_dns, V_dns, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        for n in range(nrep):
            ax0.plot(R * np.cos(th), 0.5 + n + R * np.sin(th), 'k-', lw=0.8)

        ax0.set_xlim(-L, L); ax0.set_ylim(0, nrep)
        ax0.set_aspect("equal"); ax0.set_xticks([]); ax0.set_yticks([])
        ax0.set_title(f"DNS  $t = {target_t*1e6:.0f}$ µs")
        if j == 0:
            ax0.set_ylabel("(a) DNS")

        # ── Row 1: EI homogenized ─────────────────────────────────────────────
        UL = U_hom[:nL]; UR = U_hom[nL:]
        VL = np.concatenate([UL for _ in range(nrep)])
        VR = np.concatenate([UR for _ in range(nrep)])
        ax1.tricontourf(tri_hL, VL, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        ax1.tricontourf(tri_hR, VR, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        ax1.axvline(-R, color='k', lw=1.2, ls='--')
        ax1.axvline( R, color='k', lw=1.2, ls='--')

        ax1.set_xlim(-L, L); ax1.set_ylim(0, nrep)
        ax1.set_aspect("equal"); ax1.set_xticks([]); ax1.set_yticks([])
        if j == 1:
            ax1.set_title(f"Enlarged interface  $t = {target_t*1e6:.0f}$ µs")
        else:
            ax1.set_title(f"$t = {target_t*1e6:.0f}$ µs")
        if j == 0:
            ax1.set_ylabel("(b) EI homogenized")

        # ── Row 2: 1D cross-section at y=H/2 ─────────────────────────────────
        u1d_dns = U_dns[idx_dns]
        u1d_hL  = UL[idxL]; u1d_hR = UR[idxR]

        ax2.plot(x_dns, u1d_dns, '-',  color='gray', lw=2.5, label="DNS")
        ax2.plot(xL,    u1d_hL,  'b--', lw=1.5,       label="Homogenized")
        ax2.plot(xR,    u1d_hR,  'b--', lw=1.5)
        # dotted line across the Ventcel gap [-R, R]
        ax2.plot([-R, R], [u1d_hL[-1], u1d_hR[0]], 'b:', lw=1.5,
                 label=r"Ventcel jump $[-R,R]$")

        ax2.axvline(-R, color='k', ls=':', lw=1.0)
        ax2.axvline( R, color='k', ls=':', lw=1.0)

        ax2.set_xlim(-L, L); ax2.set_ylim(-0.8, 1.1)
        ax2.set_xlabel("$x$"); ax2.grid(alpha=0.3)
        if j == 0:
            ax2.set_ylabel("Amplitude")
            ax2.legend(loc="upper left", frameon=False, fontsize=8)

    out = os.path.join(FIG, "fig_scattering_3times.pdf")
    plt.savefig(out, bbox_inches='tight')
    plt.savefig(out.replace(".pdf", ".png"), bbox_inches='tight', dpi=150)
    plt.close('all')
    print(f"  saved: {out}")
    print(f"  mu_i/mu_m={MU_I/MU_M:.1f}  rho_i/rho_m={RHO_I/RHO_M:.2f}  R={R:.3f}")
    print(f"  B={B:.3f}  S={S:.3f}  e=2R={2*R:.3f}")


if __name__ == '__main__':
    make_figure()
