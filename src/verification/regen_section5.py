#!/usr/bin/env python3
"""
Regenerate the REAL Section-5 artifacts from the verified broken-Gamma
Nitsche scheme:
  (1) convergence Table 1 + Fig 2  (MMS error vs H, real O(H^2)/O(H));
  (2) performance Table 4          (measured wall-clock + dofs + speedup).
Outputs CSV + PDF figures + LaTeX-ready rows. No hardcoded results.
"""
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _os.path.join(_SRC, 'solvers') not in _sys.path:
    _sys.path.insert(0, _os.path.join(_SRC, 'solvers'))
del _sys, _os, _SRC

import os, time, csv
import numpy as np
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from broken_nitsche_2d import assemble, MU, PI
import numpy as _np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
DATA = os.path.join(HERE, "..", "..", "data")
os.makedirs(DATA, exist_ok=True)


def solve_mms(N, timed=False):
    t0 = time.perf_counter()
    K, M, xy, bnd = assemble(N)
    x, y = xy[:, 0], xy[:, 1]
    f = MU * 2 * PI**2 * np.sin(PI * x) * np.sin(PI * y)
    b = M @ f
    free = np.setdiff1d(np.arange(K.shape[0]), bnd)
    u = np.zeros(K.shape[0])
    u[free] = spla.spsolve(K[free][:, free].tocsc(), b[free])
    dt = time.perf_counter() - t0
    ue = np.sin(PI * x) * np.sin(PI * y)
    e = u - ue
    l2 = float(np.sqrt(e @ (M @ e)))
    h1 = float(np.sqrt(e @ (K @ e) / MU))
    return dict(h=1.0 / N, ndof=K.shape[0], l2=l2, h1=h1, t=dt, nnz=K.nnz)


def rate(prev, cur, key):
    return _np.log(prev[key] / cur[key]) / _np.log(prev["h"] / cur["h"])


def main():
    print("=" * 70)
    print("  Regenerating Section 5 (real, verified pipeline)")
    print("=" * 70)
    Ns = [4, 8, 16, 32, 64]
    R = [solve_mms(N, timed=True) for N in Ns]

    # ---------- (1) convergence ----------
    print("\n  CONVERGENCE (Table 1 / Fig 2)")
    print("   H          dofs     L2 err       L2 rate   H1 err       H1 rate")
    for i, r in enumerate(R):
        rl2 = "-" if i == 0 else f"{rate(R[i-1], r, 'l2'):.2f}"
        rh1 = "-" if i == 0 else f"{rate(R[i-1], r, 'h1'):.2f}"
        print(f"  {r['h']:.4e} {r['ndof']:>7d}  {r['l2']:.4e}   {rl2:>5}   {r['h1']:.4e}   {rh1:>5}")
    hs = np.array([r["h"] for r in R])
    pL2 = np.polyfit(np.log(hs), np.log([r["l2"] for r in R]), 1)[0]
    pH1 = np.polyfit(np.log(hs), np.log([r["h1"] for r in R]), 1)[0]
    print(f"  LSQ orders: L2={pL2:.2f}  H1={pH1:.2f}")

    with open(os.path.join(DATA, "convergence_real.csv"), "w", newline="") as f:
        wc = csv.writer(f); wc.writerow(["H", "ndof", "L2", "H1"])
        for r in R:
            wc.writerow([r["h"], r["ndof"], r["l2"], r["h1"]])

    fig, ax = plt.subplots(figsize=(4.3, 3.4))
    ax.loglog(hs, [r["h1"] for r in R], "o-", label=r"$\|u-u_H\|_{H^1}$ (rate %.2f)" % pH1)
    ax.loglog(hs, [r["l2"] for r in R], "s-", label=r"$\|u-u_H\|_{L^2}$ (rate %.2f)" % pL2)
    ax.loglog(hs, [r["h1"] for r in R][0] * (hs / hs[0]), "k--", lw=0.8, label=r"$\mathcal{O}(H)$")
    ax.loglog(hs, [r["l2"] for r in R][0] * (hs / hs[0]) ** 2, "k:", lw=0.8, label=r"$\mathcal{O}(H^2)$")
    ax.set_xlabel("mesh size $H$"); ax.set_ylabel("error"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_convergence.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  -> figures/fig2_convergence.pdf (overwritten with REAL data)")

    # ---------- (2) performance (measured) ----------
    print("\n  PERFORMANCE (Table 4): measured wall-clock vs dofs")
    fine = R[-1]                       # finest = reference
    print("   H          dofs     L2 err       CPU time (s)   speedup vs finest")
    perf_rows = []
    for r in R[1:]:                    # skip coarsest
        sp = fine["t"] / r["t"]
        perf_rows.append((r["h"], r["ndof"], r["l2"], r["t"], sp))
        tag = "(reference)" if r is fine else f"{sp:.1f}x"
        print(f"  {r['h']:.4e} {r['ndof']:>7d}  {r['l2']:.4e}   {r['t']:.3f}        {tag}")

    print("\n  LaTeX rows (convergence):")
    for i, r in enumerate(R):
        rl2 = "---" if i == 0 else f"{rate(R[i-1], r, 'l2'):.2f}"
        rh1 = "---" if i == 0 else f"{rate(R[i-1], r, 'h1'):.2f}"
        print(f"    {r['h']:.4f} & {r['ndof']} & {r['h1']:.2e} & {rh1} & {r['l2']:.2e} & {rl2} \\\\")
    print("\n  done.")


if __name__ == "__main__":
    main()
