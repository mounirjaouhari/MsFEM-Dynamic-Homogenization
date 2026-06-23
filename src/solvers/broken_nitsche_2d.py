#!/usr/bin/env python3
"""
Symmetric Nitsche coupling at the SINGLE interface Gamma = {x=0.5}, on a
space that is continuous P1 within each subdomain and BROKEN only at Gamma
(Gamma nodes duplicated for the right side). This is the paper's setting and
it avoids the full-DG assembly artifact found earlier.

a_h(u,v) = sum_T mu (grad u, grad v)_T
         - int_Gamma {mu d_n u}[v] - int_Gamma [u]{mu d_n v}
         + (gamma/h) int_Gamma [u][v],
with n the +x normal, [.]=()_R-()_L, {.}=average.

Hand-rolled structured P1 assembly (full control, no black box). Verified by
MMS: u_exact=sin(pi x)sin(pi y) (continuous, =0 on dOmega, [u]=0 on Gamma);
the scheme must recover it at O(H) in H1, O(H^2) in L2.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

MU = 1.0
PI = np.pi


def build(N):
    """Structured triangulation of [0,1]^2, broken at x=1/2 (N even)."""
    assert N % 2 == 0
    nid = lambda i, j: j * (N + 1) + i
    base = (N + 1) ** 2
    gdup = {j: base + j for j in range(N + 1)}          # dup of Gamma node (N/2,j)
    ndof = base + (N + 1)
    xy = np.zeros((ndof, 2))
    for j in range(N + 1):
        for i in range(N + 1):
            xy[nid(i, j)] = (i / N, j / N)
        xy[gdup[j]] = (0.5, j / N)                       # duplicate coords

    def dofs_of(i, j_tri_nodes, right):
        out = []
        for (ii, jj) in j_tri_nodes:
            if ii == N // 2 and right:
                out.append(gdup[jj])
            else:
                out.append(nid(ii, jj))
        return out

    tris = []                                            # (global dof triple, node coords)
    for j in range(N):
        for i in range(N):
            right = i >= N // 2
            t1 = [(i, j), (i + 1, j), (i + 1, j + 1)]
            t2 = [(i, j), (i + 1, j + 1), (i, j + 1)]
            for t in (t1, t2):
                d = dofs_of(i, t, right)
                xs = np.array([(ii / N, jj / N) for (ii, jj) in t])
                tris.append((d, xs))

    # Gamma edges: between (N/2,j) and (N/2,j+1); left tri = cell(N/2-1).T1,
    # right tri = cell(N/2).T2.  Record (left dofs+coords, right dofs+coords,
    # left edge-local idx, right edge-local idx).
    gam = []
    iL, iR = N // 2 - 1, N // 2
    for j in range(N):
        # left cell T1 = [(iL,j),(iL+1,j),(iL+1,j+1)] ; edge nodes = (N/2,j),(N/2,j+1) = local 1,2
        tL = [(iL, j), (iL + 1, j), (iL + 1, j + 1)]
        dL = dofs_of(iL, tL, right=False)
        xL = np.array([(ii / N, jj / N) for (ii, jj) in tL])
        eL = [1, 2]                                      # local idx of (N/2,j),(N/2,j+1) in tL
        # right cell T2 = [(iR,j),(iR+1,j+1),(iR,j+1)] ; edge nodes = (N/2,j),(N/2,j+1)=local 0,2
        tR = [(iR, j), (iR + 1, j + 1), (iR, j + 1)]
        dR = dofs_of(iR, tR, right=True)
        xR = np.array([(ii / N, jj / N) for (ii, jj) in tR])
        eR = [0, 2]
        gam.append((dL, xL, eL, dR, xR, eR))

    bnd = [nid(i, j) for j in range(N + 1) for i in range(N + 1)
           if i in (0, N) or j in (0, N)]
    return xy, tris, gam, ndof, np.array(sorted(set(bnd)))


def p1(xs):
    """P1: area, grad of the 3 shape functions (constant)."""
    (x1, y1), (x2, y2), (x3, y3) = xs
    detT = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    A = 0.5 * abs(detT)
    b = np.array([y2 - y3, y3 - y1, y1 - y2]) / detT
    c = np.array([x3 - x2, x1 - x3, x2 - x1]) / detT
    return A, b, c                                       # grad phi_k = (b[k], c[k])


def assemble(N, gamma0=10.0):
    xy, tris, gam, ndof, bnd = build(N)
    K = sp.lil_matrix((ndof, ndof))
    M = sp.lil_matrix((ndof, ndof))
    Mloc = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]) / 12.0
    for d, xs in tris:
        A, b, c = p1(xs)
        Ke = MU * A * (np.outer(b, b) + np.outer(c, c))
        Me = A * Mloc
        for a in range(3):
            for bb in range(3):
                K[d[a], d[bb]] += Ke[a, bb]
                M[d[a], d[bb]] += Me[a, bb]
    h = 1.0 / N
    # 1D edge mass over jump nodal values (J at the two edge endpoints)
    Me1d = h * np.array([[1 / 3, 1 / 6], [1 / 6, 1 / 3]])
    for (dL, xL, eL, dR, xR, eR) in gam:
        AL, bL, cL = p1(xL)
        AR, bR, cR = p1(xR)
        gxL, gxR = bL, bR                                 # d/dx of the 3 shape fns (n=+x)
        # edge dofs: lower endpoint (j) and upper (j+1)
        dof_Lj, dof_Lb = dL[eL[0]], dL[eL[1]]
        dof_Rj, dof_Rb = dR[eR[0]], dR[eR[1]]
        edge = [dof_Lj, dof_Lb, dof_Rj, dof_Rb]
        sgn = np.array([-1.0, -1.0, +1.0, +1.0])          # [v]=v_R-v_L at endpoints
        # --- penalty (gamma/h) int [u][v] : jump nodal vals J=(j,b) ---
        Bmap = np.zeros((2, 4)); Bmap[0, 0] = -1; Bmap[0, 2] = 1; Bmap[1, 1] = -1; Bmap[1, 3] = 1
        Pen = (gamma0 * MU / h) * (Bmap.T @ Me1d @ Bmap)
        for a in range(4):
            for bb in range(4):
                K[edge[a], edge[bb]] += Pen[a, bb]
        # --- consistency: -int {mu d_n u}[v] (and transpose) ---
        # {mu d_n u} constant = 0.5*mu*(sum_k u_LT gxL + sum_k u_RT gxR)
        # int_Gamma [v] = (h/2)*sum(sgn * v_edge)
        u_cols = list(dL) + list(dR)                      # 6 trial dofs
        u_coef = 0.5 * MU * np.concatenate([gxL, gxR])    # coefficient of each on {mu d_n u}
        for a in range(4):                                # v edge dof
            w_v = (h / 2.0) * sgn[a]
            for q in range(6):
                val = -w_v * u_coef[q]
                K[edge[a], u_cols[q]] += val              # -{mu d_n u}[v]
                K[u_cols[q], edge[a]] += val              # -[u]{mu d_n v} (transpose)
    return K.tocsr(), M.tocsr(), xy, bnd


def mms(N, gamma0=10.0):
    K, M, xy, bnd = assemble(N, gamma0)
    x, y = xy[:, 0], xy[:, 1]
    f = MU * 2 * PI**2 * np.sin(PI * x) * np.sin(PI * y)
    b = M @ f
    free = np.setdiff1d(np.arange(K.shape[0]), bnd)
    u = np.zeros(K.shape[0])
    u[free] = spla.spsolve(K[free][:, free].tocsc(), b[free])
    ue = np.sin(PI * x) * np.sin(PI * y)
    e = u - ue
    l2 = np.sqrt(e @ (M @ e))                              # mass-weighted L2
    h1 = np.sqrt(e @ (K @ e) / MU)                         # energy seminorm ~ H1
    return 1.0 / N, l2, h1


def main():
    print("=" * 64)
    print("  Broken-at-Gamma Nitsche : MMS convergence (gamma0=10)")
    print("=" * 64)
    rows = [mms(N) for N in (4, 8, 16, 32)]
    hs = np.array([r[0] for r in rows]); el2 = np.array([r[1] for r in rows]); eh1 = np.array([r[2] for r in rows])
    print("    h          L2 err       L2 rate    H1 err       H1 rate")
    for i in range(len(rows)):
        rl2 = "" if i == 0 else f"{np.log(el2[i-1]/el2[i])/np.log(hs[i-1]/hs[i]):.2f}"
        rh1 = "" if i == 0 else f"{np.log(eh1[i-1]/eh1[i])/np.log(hs[i-1]/hs[i]):.2f}"
        print(f"  {hs[i]:.4e}  {el2[i]:.4e}  {rl2:>6}    {eh1[i]:.4e}  {rh1:>6}")
    pL2 = np.polyfit(np.log(hs), np.log(el2), 1)[0]
    pH1 = np.polyfit(np.log(hs), np.log(eh1), 1)[0]
    print(f"  LSQ orders: L2={pL2:.3f} (expect ~2), H1={pH1:.3f} (expect ~1)")
    ok = (1.8 <= pL2 <= 2.2) and (0.85 <= pH1 <= 1.2)
    print("-" * 64)
    print(f"  broken-Gamma Nitsche MMS: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
