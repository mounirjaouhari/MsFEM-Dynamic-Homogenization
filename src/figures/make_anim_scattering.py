#!/usr/bin/env python3
"""
Animated transient scattering: DNS vs Enlarged-Interface (anim_scattering.gif).

Three panels animated over time:
  Top    — DNS 2D wavefield (circular inclusion, FEM resolved)
  Middle — EI homogenized (two half-domains + Ventcel jump conditions)
  Bottom — 1D cross-section at y=H/2 (DNS vs EI + Ventcel gap)

Same physics as make_fig_scattering.py (structural regime):
  mu_i/mu_m=6.5, rho_i/rho_m=3.1, R=0.282, B=-0.31, S=0.53, e=2R

Output: figures/anim_scattering.gif  (~65 frames, 15 fps)
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
import matplotlib.animation as animation
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm
from skfem.helpers import dot, grad
from interface_params_2d import MU_M, MU_I, RHO_M, RHO_I, R

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)

L        = 3.0
H_period = 1.0
C_M      = np.sqrt(MU_M / RHO_M)
B        = -0.31
S        =  0.53
N_FRAMES =  65
FPS      =  15

def incl(x, y):
    return x**2 + (y - 0.5)**2 <= R**2

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
    return m, stiff.assemble(vb), mass.assemble(vb)

def get_homog_matrices():
    e_gap = 2.0 * R
    nx, ny = 120, 31
    mL = MeshTri.init_tensor(np.linspace(-L, -e_gap/2, nx//2), np.linspace(0.0, 1.0, ny))
    mR = MeshTri.init_tensor(np.linspace( e_gap/2,  L, nx//2), np.linspace(0.0, 1.0, ny))
    vbL = Basis(mL, ElementTriP1()); vbR = Basis(mR, ElementTriP1())
    @BilinearForm
    def stiff_hom(u, v, w): return MU_M * dot(grad(u), grad(v))
    @BilinearForm
    def mass_hom(u, v, w):  return RHO_M * u * v
    KL = stiff_hom.assemble(vbL); ML = mass_hom.assemble(vbL)
    KR = stiff_hom.assemble(vbR); MR = mass_hom.assemble(vbR)
    nL = KL.shape[0]
    def _iface(mesh, x0):
        p = mesh.p.T; n = np.where(np.abs(p[:,0]-x0)<1e-9)[0]
        o = np.argsort(p[n,1]); return n[o], np.gradient(p[n[o],1])
    rL, dL = _iface(mL, -e_gap/2); rR, dR = _iface(mR, e_gap/2)
    dy = 0.5*(dL+dR)
    Kh = sp.bmat([[KL,None],[None,KR]], format="lil")
    Mh = sp.bmat([[ML,None],[None,MR]], format="lil")
    B_enl = B + e_gap/H_period; S_enl = S + e_gap/H_period
    kap = MU_M/(H_period*B_enl); sm = H_period*S_enl*RHO_M
    for a, b, d in zip(rL, rR+nL, dy):
        Kh[a,a]+=kap*d; Kh[b,b]+=kap*d; Kh[a,b]-=kap*d; Kh[b,a]-=kap*d
        Mh[a,a]+=0.25*sm*d; Mh[b,b]+=0.25*sm*d
        Mh[a,b]+=0.25*sm*d; Mh[b,a]+=0.25*sm*d
    class Hm:
        def __init__(self,a,b): self.mL=a; self.mR=b
    return Hm(mL,mR), Kh.tocsc(), Mh.tocsc(), nL

def _leapfrog_frames(K, M, x_coords, save_times):
    ml = np.asarray(M.sum(axis=1)).ravel(); inv_ml = 1.0/ml
    lam = spla.eigsh(K.tocsc(),k=1,M=sp.diags(ml),which="LA",return_eigenvectors=False)[0]
    dt  = 0.9*2.0/np.sqrt(lam)
    x0, sig = -1.5, 0.25
    U0 = np.exp(-(x_coords-x0)**2/(2*sig**2))
    V0 = C_M*(x_coords-x0)/(sig**2)*U0
    acc0 = inv_ml*(-(K@U0))
    U_curr = U0.copy(); U_prev = U0 - dt*V0 + 0.5*dt**2*acc0
    frames = []; t = 0.0
    for target_t in save_times:
        while t < target_t:
            acc = inv_ml*(-(K@U_curr))
            U_next = 2*U_curr - U_prev + dt**2*acc
            U_prev, U_curr = U_curr, U_next
            t += dt
        frames.append(U_curr.copy())
    return frames

def _tiled_tri(mesh, nrep=3):
    p = mesh.p; bt = mesh.t.T
    xs,ys,ts = [],[],[]
    for n in range(nrep):
        xs.append(p[0]); ys.append(p[1]+n); ts.append(bt+n*p.shape[1])
    return mtri.Triangulation(np.concatenate(xs),np.concatenate(ys),np.concatenate(ts,axis=0))

def make_animation():
    print("  Assembling DNS matrices...")
    m_dns, K_dns, M_dns = get_dns_matrices()
    print("  Assembling homogenized matrices (Ventcel jump)...")
    m_hom, K_hom, M_hom, nL = get_homog_matrices()

    t_hit   = abs(-1.5) / C_M
    T_total = 2.5 * t_hit
    save_times = np.linspace(0.0, T_total, N_FRAMES)

    p_hom_x = np.concatenate([m_hom.mL.p[0], m_hom.mR.p[0]])

    print("  Running DNS lockstep...")
    frames_dns = _leapfrog_frames(K_dns, M_dns, m_dns.p[0], save_times)
    print("  Running EI homogenized lockstep...")
    frames_hom = _leapfrog_frames(K_hom, M_hom, p_hom_x, save_times)

    # triangulations
    nrep    = 3
    tri_dns = _tiled_tri(m_dns, nrep)
    tri_hL  = _tiled_tri(m_hom.mL, nrep)
    tri_hR  = _tiled_tri(m_hom.mR, nrep)

    # 1D slice nodes at y=0.5
    p_dns = m_dns.p.T
    idx_d = np.where(np.abs(p_dns[:,1]-0.5)<1e-6)[0]
    idx_d = idx_d[np.argsort(p_dns[idx_d,0])]; x_d = p_dns[idx_d,0]
    pL = m_hom.mL.p.T; pR = m_hom.mR.p.T
    idxL = np.where(np.abs(pL[:,1]-0.5)<1e-6)[0]
    idxL = idxL[np.argsort(pL[idxL,0])]; xL = pL[idxL,0]
    idxR = np.where(np.abs(pR[:,1]-0.5)<1e-6)[0]
    idxR = idxR[np.argsort(pR[idxR,0])]; xR = pR[idxR,0]

    lim = 0.4; lvl = np.linspace(-lim, lim, 41)
    th  = np.linspace(0.0, 2*np.pi, 60)

    fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.5),
                              gridspec_kw={'height_ratios': [1, 1, 1.2]})

    def animate(i):
        ax0, ax1, ax2 = axes
        ax0.clear(); ax1.clear(); ax2.clear()

        U_dns = frames_dns[i]
        U_hom = frames_hom[i]
        t_us  = save_times[i] * 1e6

        # DNS
        V_dns = np.concatenate([U_dns for _ in range(nrep)])
        ax0.tricontourf(tri_dns, V_dns, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        for n in range(nrep):
            ax0.plot(R*np.cos(th), 0.5+n+R*np.sin(th), 'k-', lw=0.8)
        ax0.set_xlim(-L, L); ax0.set_ylim(0, nrep)
        ax0.set_aspect("equal"); ax0.set_xticks([]); ax0.set_yticks([])
        ax0.set_title(f"DNS  —  $t = {t_us:.1f}$ µs", fontsize=9)

        # EI homogenized
        UL = U_hom[:nL]; UR = U_hom[nL:]
        VL = np.concatenate([UL for _ in range(nrep)])
        VR = np.concatenate([UR for _ in range(nrep)])
        ax1.tricontourf(tri_hL, VL, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        ax1.tricontourf(tri_hR, VR, levels=lvl, cmap=pubstyle.FIELD_CMAP, extend="both")
        ax1.axvline(-R, color='k', lw=1.2, ls='--')
        ax1.axvline( R, color='k', lw=1.2, ls='--')
        ax1.set_xlim(-L, L); ax1.set_ylim(0, nrep)
        ax1.set_aspect("equal"); ax1.set_xticks([]); ax1.set_yticks([])
        ax1.set_title(f"Enlarged interface  ($e=2R$, $B={B}$, $S={S}$)  —  $t = {t_us:.1f}$ µs", fontsize=9)

        # 1D cross-section
        u1d_d  = U_dns[idx_d]
        u1d_hL = UL[idxL]; u1d_hR = UR[idxR]
        ax2.plot(x_d, u1d_d,  '-',  color='gray', lw=2.5, label="DNS")
        ax2.plot(xL,  u1d_hL, 'b--', lw=1.5,       label="Homogenized")
        ax2.plot(xR,  u1d_hR, 'b--', lw=1.5)
        ax2.plot([-R, R], [u1d_hL[-1], u1d_hR[0]], 'b:', lw=1.5, label=r"Ventcel jump $[-R,R]$")
        ax2.axvline(-R, color='k', ls=':', lw=1.0)
        ax2.axvline( R, color='k', ls=':', lw=1.0)
        ax2.set_xlim(-L, L); ax2.set_ylim(-1.1, 1.1)
        ax2.set_xlabel("$x$"); ax2.set_ylabel(r"Amplitude $U(x,\,y=0.5)$")
        ax2.legend(loc="upper right", frameon=False, fontsize=8)
        ax2.grid(alpha=0.3)

        fig.tight_layout()

    print(f"  Rendering {N_FRAMES} frames...")
    ani = animation.FuncAnimation(fig, animate, frames=N_FRAMES, blit=False)
    out = os.path.join(FIG, "anim_scattering.gif")
    ani.save(out, writer=animation.PillowWriter(fps=FPS), dpi=80)
    plt.close(fig)
    print(f"  saved: {out}")

if __name__ == '__main__':
    make_animation()
