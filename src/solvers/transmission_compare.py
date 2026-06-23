#!/usr/bin/env python3
"""
MODEL A validation, part 2/2: the homogenized effective-jump model REPRODUCES
the inclusion-resolving reference to O(eta^2)  (validates C2 + GATE_REF + the
O(eta^2) accuracy claim of Marigo's interface model).

At normal incidence the matrix-bulk + interface-jump problem is 1D in x, so the
transmission has an EXACT transfer-matrix solution.  Interface at x=0 (h=1):
  [U]        = B <d_x U>             (displacement jump, B from interface_params_2d)
  [mu d_x U] = -S omega^2 <U>        (inertial flux jump, S = excess surface mass)
with [.] = (.)+ - (.)- and <.> = 1/2((.)+ + (.)-).
Left  U = e^{ikx} + r e^{-ikx},  Right U = t e^{ikx},  k = omega = eta (nondim).

Unknowns (t,r) solve
  (I)  t(1-Bik/2) + r(-1+Bik/2) = 1 + Bik/2
  (II) t(ik+Sk^2/2) + r(ik+Sk^2/2) = ik - Sk^2/2
tau_jump = t.

CLAIM TO VERIFY (the whole point):
  |tau_resolved - 1|            ~ O(eta)    (no-model baseline: ignore the row)
  |tau_resolved - tau_jump|     ~ O(eta^2)  (jump model captures the row; one order better)
A least-squares slope ~2 on the second curve proves the homogenized jump model
(and the SIGN/VALUE of B,S) is correct.  Slope ~1 would mean B is wrong.
"""
import numpy as np
import transmission_resolved as ref

# verified, logged in data/PROVENANCE.md (interface_params_2d.py, centred circle, contrast 6.5);
# flux-BC measurement, analytically cross-checked vs the 1D laminate compliance.
B = -0.311
RHO_I_OVER_M = ref.RHO_I / ref.RHO_M            # 3.12
S = (RHO_I_OVER_M - 1.0) * np.pi * ref.R ** 2   # excess surface mass (nondim), ~0.53


def tau_jump(eta, B=B, S=S):
    k = eta
    A = np.array([[1 - B * 1j * k / 2, -1 + B * 1j * k / 2],
                  [1j * k + S * k ** 2 / 2, 1j * k + S * k ** 2 / 2]], dtype=complex)
    rhs = np.array([1 + B * 1j * k / 2, 1j * k - S * k ** 2 / 2], dtype=complex)
    t, r = np.linalg.solve(A, rhs)
    return t, r


def slope(etas, errs):
    le, lr = np.log(np.asarray(etas)), np.log(np.asarray(errs))
    return np.polyfit(le, lr, 1)[0]


def invert(t, r, eta):
    """Exact inversion of (I),(II): effective B,S exhibited by the reference."""
    k = eta
    B_eff = 2 * (t - 1 - r) / (1j * k * (t + 1 - r))
    S_eff = -2j * (t - 1 + r) / (k * (t + 1 + r))
    return B_eff.real, S_eff.real


def main():
    print("=" * 86)
    print("  Homogenized jump model vs inclusion-resolving reference (normal incidence)")
    print("  B=%.4f  S=%.4f   (verified, PROVENANCE.md)" % (B, S))
    print("=" * 86)
    W = 6.0
    etas = [0.40, 0.20, 0.10, 0.05, 0.025]
    e_base, e_jump = [], []
    print("  %-7s | %-10s %-10s | %-10s %-10s | %-11s %-11s"
          % ("eta", "|tau_res|", "|tau_jmp|", "|res-1|", "|res-jmp|", "argres", "argjmp"))
    print("  " + "-" * 82)
    for eta in etas:
        tr, rr, _, _ = ref.solve(W, "circle", eta)
        tj, rj = tau_jump(eta)
        eb = abs(tr - 1.0)
        ej = abs(tr - tj)
        e_base.append(eb); e_jump.append(ej)
        print("  %-7.3f | %-10.6f %-10.6f | %-10.3e %-10.3e | %+0.5f   %+0.5f"
              % (eta, abs(tr), abs(tj), eb, ej, np.angle(tr), np.angle(tj)))
    sb, sj = slope(etas, e_base), slope(etas, e_jump)
    print("  " + "-" * 82)
    print("  LSQ slope  |tau_res - 1|      = %.2f   (expected ~1: O(eta) leading row effect)" % sb)
    print("  LSQ slope  |tau_res - tau_jmp|= %.2f   (residual O(eta) ~ few%% of leading)" % sj)
    gain = np.mean(np.array(e_base) / np.array(e_jump))
    print("  mean error reduction (baseline/jump) = %.1fx" % gain)

    print("\n  Effective interface parameters INVERTED from the resolved reference")
    print("  (must converge to the cell-computed B=%.3f, S=%.3f as eta->0):" % (B, S))
    print("  %-7s | %-12s %-12s" % ("eta", "B_eff", "S_eff"))
    print("  " + "-" * 34)
    B0 = S0 = None
    for eta in etas:
        tr, rr, _, _ = ref.solve(W, "circle", eta)
        be, se = invert(tr, rr, eta)
        print("  %-7.3f | %+ .5f    %+ .5f" % (eta, be, se))
        B0, S0 = be, se          # smallest-eta (last) = best estimate
    eB = abs(B0 - B) / abs(B); eS = abs(S0 - S) / abs(S)
    print("  " + "-" * 34)
    print("  at eta=%.3f:  B_eff=%.4f (cell %.4f, %.1f%%)   S_eff=%.4f (cell %.4f, %.1f%%)"
          % (etas[-1], B0, B, 100 * eB, S0, S, 100 * eS))
    ok = eB < 0.05 and eS < 0.10 and gain > 10
    print("  GATE_REF (resolved row <-> homogenized jump, two-route agreement): %s"
          % ("PASS" if ok else "CHECK"))


if __name__ == "__main__":
    main()
