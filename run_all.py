#!/usr/bin/env python3
"""
Master orchestration script — MsFEM-Dynamic-Homogenization
Runs all verification gates and figure scripts in the correct order.

Usage:
    python run_all.py            # full run (all gates + all figures)
    python run_all.py --gates    # gates only (no figures, ~2 min)

Expected total runtime: 5–15 minutes depending on hardware.
All figures are written to ../figures/.
"""
import subprocess
import sys
import time
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, label):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  python {script}")
    print(f"{'='*60}")
    t0 = time.perf_counter()
    result = subprocess.run(
        [sys.executable, os.path.join(HERE, script)],
        capture_output=False,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    status = "PASS" if result.returncode == 0 else "FAIL"
    print(f"\n  [{status}] {label} ({elapsed:.1f}s)")
    return status, elapsed


def main():
    gates_only = "--gates" in sys.argv

    results = []

    # ------------------------------------------------------------------
    # PHASE 1 — Verification gates (must all PASS before figures)
    # ------------------------------------------------------------------
    print("\n" + "#"*60)
    print("  PHASE 1 — Verification gates")
    print("#"*60)

    gates = [
        (os.path.join("src", "verification", "smoke_poisson_mms.py"),     "GATE_SMOKE : Poisson MMS O(H^2)/O(H)"),
        (os.path.join("src", "solvers",      "cell_problem_2d.py"),        "GATE_VNV[cell] : mu* within Hashin-Shtrikman"),
        (os.path.join("src", "verification", "nitsche_coercivity_2d.py"), "GATE_VNV[Nitsche] : symmetry + coercivity"),
        (os.path.join("src", "solvers",      "transient_energy_2d.py"),   "GATE_VNV[time] : energy drift < 1e-14"),
        (os.path.join("src", "solvers",      "broken_nitsche_2d.py"),     "GATE_VNV[MMS] : broken-Gamma Nitsche O(H^2)/O(H)"),
        (os.path.join("src", "verification", "regen_section5.py"),        "GATE_COMPUTE : Table 1 + Table 4"),
    ]
    for script, label in gates:
        status, elapsed = run(script, label)
        results.append((label, status, elapsed))

    if gates_only:
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # PHASE 2 — Interface parameters (Table 2)
    # ------------------------------------------------------------------
    print("\n" + "#"*60)
    print("  PHASE 2 — Interface parameters")
    print("#"*60)

    status, elapsed = run(os.path.join("src", "solvers", "interface_params_2d.py"),
                         "Interface parameters B, B2, C, C1, S (Table 2)")
    results.append(("Interface parameters", status, elapsed))

    # ------------------------------------------------------------------
    # PHASE 3 — Figures (depend on Phase 1 and 2)
    # ------------------------------------------------------------------
    print("\n" + "#"*60)
    print("  PHASE 3 — Figures")
    print("#"*60)

    figures = [
        (os.path.join("src", "figures", "make_meshfig.py"),     "Fig 2 — mesh hierarchy"),
        (os.path.join("src", "figures", "make_basisfig.py"),    "Fig 3 — spectral MsFEM basis functions"),
        (os.path.join("src", "figures", "make_cellfig.py"),     "Fig 4 — strip-cell corrector fields"),
        (os.path.join("src", "figures", "make_figs_modelA.py"), "Fig 5 — transmission validation (gate/ref)"),
        (os.path.join("src", "figures", "make_panel6.py"),      "Fig 6 — 6-panel solution"),
        (os.path.join("src", "figures", "make_fig_energy.py"),  "Fig 8 — discrete energy conservation"),
        (os.path.join("src", "solvers", "resonant_mass_2d.py"), "Fig 9 — resonant band gap"),
        (os.path.join("src", "figures", "sensitivity_wave.py"), "Sensitivity + wave snapshot figures"),
    ]
    for script, label in figures:
        status, elapsed = run(script, label)
        results.append((label, status, elapsed))

    _print_summary(results)


def _print_summary(results):
    print("\n" + "="*60)
    print("  SUMMARY")
    print("="*60)
    n_pass = sum(1 for _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    for label, status, elapsed in results:
        mark = "+" if status == "PASS" else "x"
        print(f"  {mark} {status:4s}  ({elapsed:5.1f}s)  {label}")
    print(f"\n  Total: {n_pass} PASS, {n_fail} FAIL")
    if n_fail:
        print("  Some steps FAILED — check output above.")
        sys.exit(1)
    else:
        print("  All steps PASSED.")


if __name__ == "__main__":
    main()
