# PROVENANCE LEDGER — MsFEM 2D shear-wave paper (real S2 computation)

Every numeric claim destined for the manuscript is logged here with the command that
produced it. No row ⇒ the claim is deleted (GATE 0). Stack: scikit-fem 12.0.1 (FEniCSx
unavailable on this Windows host; SWITCH_FEM_STACK per GATE_SMOKE rule). Python 3.13.

## Gate ledger
- **GATE_SMOKE**: PASS — C-smoke-*.
- **GATE_VNV**: PASS — cell (μ* in HS), time (energy 5e-15), Nitsche (broken-Γ MMS O(H²)/O(H)).
- **GATE_COMPUTE**: PASS — every §5 number (Tables 1,4; Figs 2,3,4,5; μ*,ρ*) is produced by this pipeline and logged here; no fabricated value remains in main.tex.
- **GATE_REF**: PARTIAL — verification uses a manufactured solution (MMS) as reference, which is rigorous for the scheme; an *inclusion-resolving* fine reference (numerical Term-I homogenization error) is deferred to a companion study (stated in §5).
- **GATE_NOGO**: CONFIRM/REFRAME — honest results support the corrected contributions (spectral enrichment; Nitsche jumps; O(η+H); exact energy conservation; verified numerics). The microstructure-capturing speedup vs a resolved reference is scoped as future work.

## §5 regeneration (real, replaces fabricated)
| claim_id | value | command | artifact |
|---|---|---|---|
| C-conv-table | H1 rate→1.01, L2 rate→1.99 (5 levels) | `python regen_section5.py` | data/convergence_real.csv, figures/fig2_convergence.pdf, Table 1 |
| C-perf | measured s; speedup 3.8×–60× vs finest | `python regen_section5.py` | Table 4 |
| C-sensitivity | μ*/μ_m: 1.00→1.66 over contrast 1→40, all within HS | `python sensitivity_wave.py` | figures/fig4_sensitivity.pdf |
| C-wave | real field snapshot U(x,1/2), max|U|=0.16 | `python sensitivity_wave.py` | figures/fig5_wave_profiles.pdf |

## Claims

| claim_id | value | metric definition | command | artifact | maps_to_method |
|---|---|---|---|---|---|
| C-smoke-L2-order | 1.968 | LSQ slope of log(L2 error) vs log(h), P1, MMS u=sin(πx)sin(πy) | `python smoke_poisson_mms.py` | stdout | GATE_SMOKE: real 2D scikit-fem Poisson solve |
| C-smoke-H1-order | 0.988 | LSQ slope of log(H1 error) vs log(h), same | `python smoke_poisson_mms.py` | stdout | GATE_SMOKE |
| C-smoke-deterministic | true | identical L2 errors on repeated run | `python smoke_poisson_mms.py` | stdout | GATE_SMOKE |
| C-mustar | 17.47 GPa (iso) | μ*_11=μ*_22=1.747e10, μ*_12≈-1.7e7 (≈0), SPD, 2D cell problem h=1/128 | `python cell_problem_2d.py` | stdout | GATE_VNV[cell]: real periodic cell solve |
| C-mustar-bounds | within HS | Reuss 15.22 < HS- 17.39 ≤ μ* 17.47 ≤ HS+ 22.64 < Voigt 28.5 GPa (matrix-continuous → near HS-) | `python cell_problem_2d.py` | stdout | GATE_VNV[cell] |
| C-rhostar | 3825 kg/m³ | ρ* = (1-vf)ρ_m + vf ρ_i = 0.75·2500 + 0.25·7800 | analytic (matches ⟨ρ⟩_Y) | GATE_VNV[cell] |
| C-energy-drift | 5.1e-15 | max_t \|E_H(t)/E_H(0) − 1\|, leapfrog modified energy, 4000 steps, f=0 | `python transient_energy_2d.py` | GATE_VNV[time]: REAL solve (replaces fabricated np.ones) |
| C-energy-genuine | drift>0 | roundoff-level, not exactly constant ⇒ not fabricated | `python transient_energy_2d.py` | GATE_VNV[time] |
| FIG-energy | figures/fig3_energy_real.pdf | E_H(t)/E_H(0) + log drift, real 4000-step solve | `python make_fig_energy.py` | candidate replacement for fabricated Fig.3 |
| C-nitsche-symmetric | 0.0 | \|K_H−K_Hᵀ\|/\|K_H\| for symmetric interior-penalty form | `python nitsche_coercivity_2d.py` | GATE_VNV[Nitsche]: exact symmetry — VERIFIED |
| C-nitsche-gammacrit | ≈2.0 (refine=2) | at refine=2: SPD with λ_min=0.42, indefinite below γ₀≈2 (inverse-trace threshold) | `python nitsche_coercivity_2d.py` (refine=2) | GATE_VNV[Nitsche]: coercivity threshold VERIFIED at coarse level |

> **GATE_VNV[Nitsche] = PARTIAL / OPEN (honest status).** The symmetric Nitsche form is *exactly symmetric* and *coercive (SPD)* at `refine=2` with `γ_crit≈2`. BUT the assembled DG (`ElementTriDG`) operator becomes **numerically singular at refine≥3** (λ_min → ~1e-17), which made `convergence_nitsche_mms.py` return NaN. This is an unresolved scikit-fem DG facet-assembly issue, NOT a verified physical result. The earlier "PASS" was overstated and is retracted. → The interface-convergence study (GATE_REF / GATE_COMPUTE) is **blocked** until this is fixed (debug skfem DG, or use continuous-subdomain + single-interface Nitsche instead of full DG).
>
> Note: the **discretization order (Term II, O(H)/O(H²))** is independently verified by the continuous-Galerkin MMS smoke test (C-smoke-* above); the paper's MsFEM is conforming within each subdomain, so the single interface is a lower-order correction.
>
> **DIAGNOSIS (2026-06-16):** the singularity is a genuine nullmode at refine≥3, λ_min≈0 *independent of γ₀* (tested γ₀=10/100/1000 → all ~1e-17), i.e. a scikit-fem `ElementTriDG`+`InteriorFacetBasis` assembly artifact that scales with the interior-facet count, NOT a coercivity-threshold issue. refine=2 is clean (λ_min=0.42). → ABANDONED full-DG.
>
> **RESOLVED (2026-06-16):** implemented `broken_nitsche_2d.py` — continuous P1 within each subdomain, broken ONLY at Γ (Γ nodes duplicated for the right side), with hand-rolled symmetric Nitsche edge integrals on Γ. **MMS-VERIFIED at optimal order** → GATE_VNV[Nitsche] = **PASS** (genuine).

| C-nitsche-mms-L2 | rate→1.98 | broken-Γ Nitsche MMS, L2 rates 1.73/1.93/1.98 (LSQ 1.885) | `python broken_nitsche_2d.py` | GATE_VNV[Nitsche]: scheme coercive & convergent (VERIFIED) |
| C-nitsche-mms-H1 | rate→1.05 | H1 rates 1.31/1.15/1.05 (LSQ 1.171) | `python broken_nitsche_2d.py` | GATE_VNV[Nitsche] |

---

## MODEL REFRAME (2026-06-16): MODEL A — single row of elastic inclusions (Marigo-faithful)

Author decision (locked, see `MODEL_SPEC.md`): the paper models a **single row** of inclusions in a
matrix (not a bulk periodic composite). Bulk = matrix on both sides (ρ*=ρ_m, μ*=μ_m); all microstructure
physics is in the effective interface jump. The `Y=[0,1]²` periodic corrector + μ*_ij (old `cell_problem_2d.py`)
is **retired** as the model's bulk (μ* survives only as the Nitsche penalty scale). Spec reconstructed and
source-verified vs Pham-Maurel-Marigo PRSA 2021, Cornaggia-Touboul-Bellis CRAS 2022, Lombard-Maurel-Marigo,
Marigo-Maurel PRSA 2016 (open access); target MMP-S J.Elasticity 2017 abstract only (paywalled) — disclosed.

### Verified interface jump parameters (real strip cell solves)
| claim_id | value | command | check |
|---|---|---|---|
| C-iface-matrix | B=B₂=C=C₁=0 (machine zero) | `python interface_params_2d.py` | control: no inclusion ⇒ no jump ✓ |
| C-iface-1D | laminate: B₂=C=C₁=0 **exactly**, B=−0.476 (flux BC; **was −0.203 pre-fix, see C-iface-B-analytic**) | `python interface_params_2d.py` | transverse params vanish in 1D (the 2D premise) ✓ |
| C-iface-circle | centred circle: B=−0.311, C/μ_m=−0.914; B₂,C₁≈1e-3≈0 (flux BC; **B was −0.130 pre-fix, see C-iface-B-fix**) | `python interface_params_2d.py` | B,C≠0; B₂=C₁=0 by reflection symmetry ✓ |
| C-iface-B2C1 | tilted ellipse: B₂=−0.04676, C₁=+0.04676, **B₂=−C₁ to machine precision** (flux BC; superseded the stale "+0.048, 2.3%" pre-fix value; see C-iface-B2C1-exact) | `python interface_params_2d.py` | sourced energy-consistency identity (Cornaggia Rem.1) verified NONTRIVIALLY ✓ |
| C-iface-S | excess surface mass (ρ_i−ρ_m)·πR² = 1325 (dimensionless cell) | `python interface_params_2d.py` | static inertial parameter, O(h); NOT (ρ*/ρ_m)(a/h) |

> Definitions (MODEL_SPEC §6): static strip correctors `div(μ(∇Q^{(j)}+e_j))=0` on `(−L,L)×(0,1)`,
> y₂-periodic; `B,B₂` = far-field plateau jumps of `Q^{(1)},Q^{(2)}`; `C,C₁ = ∫μ∂_{y2}Q^{(2,1)}`. Stable in L.
> Replaces the WRONG `B=|Y|⁻¹∫χ¹` (≡0) and the WRONG `S=(ρ*/ρ_m)(a/h)` of the old manuscript.

### Resonant surface mass B_S(ω) + band gap (metamaterial payoff; real eigenproblem)
Locally-resonant set: soft dense inclusion μ_i=0.5 GPa, ρ_i=9000, in matrix μ_m=12 GPa, ρ_m=2500.
| claim_id | value | command | check |
|---|---|---|---|
| C-res-eig | FE ω₁ matches exact clamped-disk j₀,₁·c_i/R to **5.6e-4**; dipole pair matches j₁,₁ | `python resonant_mass_2d.py` | eigensolver verified vs analytic Bessel spectrum ✓ |
| C-res-participation | a₁=0.489 (breathing); dipole modes a_n≈1e-32 | `python resonant_mass_2d.py` | only m=0 breathing modes couple to uniform field (correct physics) ✓ |
| C-res-bandgap | M_eff(ω)<0 on [1.134e5,1.427e5] rad/s = **[18.1,22.7] kHz**, opens at f₁=18.1 kHz | `python resonant_mass_2d.py` | negative-effective-mass band gap above the fundamental resonance ✓ |

> Effective dynamic surface mass `M_eff(ω)=m_s+Σ_n a_n ω_n²/(ω_n²−ω²)` (locally-resonant interface
> closure of Pham–Maurel–Marigo 2017, the source cited in the manuscript; cf. Milton 2002 for the
> general negative-effective-mass mechanism); `ω_n,φ_n` from the clamped-inclusion eigenproblem `−div(μ∇φ)=ρω²φ, φ|_{∂D}=0`; `a_n=(∫_D ρφ_n)²`
> with mass-normalized modes. This is the genuine frequency-dependent/singular `S(ω)` the old manuscript's
> `S=(ρ*/ρ_m)(a/h)` failed to represent. (m_s static-mass normalization to be finalized at §5 integration.)

### CORRECTION: normal compliance B (flux-BC), analytically verified
The earlier Dirichlet-plateau B measurement UNDER-measured B (the perturbation is forced to 0 at the
truncation ends). Corrected to an **inhomogeneous-Neumann (matrix-flux) BC** so the corrector genuinely
plateaus; B measured at the ends.
| claim_id | value | command | check |
|---|---|---|---|
| C-iface-B-fix | B(circle)=−0.311 (was −0.130, wrong) | `python interface_params_2d.py` | flux BC, stable in L ✓ |
| C-iface-B-analytic | B(1D laminate)=−0.476 = exact 2R(μ_m/μ_i−1)=−0.477 | `python interface_params_2d.py` | analytic slab-compliance cross-check ✓ |
| C-iface-B2C1-exact | B₂=−C₁ to machine precision (consistent BCs) | `python interface_params_2d.py` | discrete energy-consistency reciprocity ✓ |

### GATE_REF: homogenized jump model ↔ inclusion-resolving reference (transmission)
Time-harmonic transmission across the row, normal incidence. Reference = resolved circular inclusions
(`transmission_resolved.py`, ABC exact for 0th order: matrix-only gives |τ|=1 to 1e-6 ✓). Homogenized
jump model = exact 1D transfer matrix with [U]=B⟨∂_xU⟩, [μ∂_xU]=−Sω²⟨U⟩ (`transmission_compare.py`).
| claim_id | value | command | check |
|---|---|---|---|
| C-trans-ref | resolved row: |τ−1|, |ρ| = O(η) (slope 1.00) | `python transmission_resolved.py` | inclusion-resolving reference (GATE_REF) ✓ |
| C-trans-model | jump model reduces τ-error vs resolved by **21×**; phase captured | `python transmission_compare.py` | homogenized jump reproduces the row ✓ |
| C-trans-tworoute | B,S inverted from reference → cell values: B_eff=−0.300 (cell −0.311, 3.5%), S_eff=0.513 (cell 0.530, 3.1%) at η=0.025 | `python transmission_compare.py` | **two-route agreement ~3%** (residual O(η)+discr.) → GATE_REF PASS |

> This closes the biggest review finding (numerics never tested the actual method/jump against a resolved
> reference). The model + its parameters are now validated against the fine-scale truth to ~3%.

### Oblique incidence (transverse content) — HONEST status
| claim_id | value | command | check |
|---|---|---|---|
| C-obl-solver | matrix-only oblique (45°, Bloch BC): |τ|=1 to 1e-6 | `python transmission_oblique.py` | Bloch/ABC oblique solver verified ✓ |
| C-obl-C-open | system-level C extraction NOT clean (subdominant) | `python transmission_oblique.py` | OPEN: C enters 0th-order transmission only at O(η²)~Ck_y²; cannot isolate from other 2nd-order terms; exact scalar stress-jump form not sourceable (paywall) |

> The 2D transverse content is established at the **cell level** (C-iface-circle: C=−0.91≠0; C-iface-B2C1-exact:
> B₂=−C₁; C-iface-1D: all transverse params vanish in 1D). A quantitative *system-level* C validation via
> oblique transmission is left as honest open work (do not overclaim in §5).

### §5 manuscript figures (Model A)
| claim_id | value | command | artifact |
|---|---|---|---|
| FIG-gateref | jump model ÷21 error vs resolved; B_eff,S_eff → cell to 3% | `python make_figs_modelA.py` | figures/fig_gateref.pdf |
| FIG-bandgap | S(ω) negative over [1,1.26]ω₁; f₁≈18 kHz | `python make_figs_modelA.py` | figures/fig_bandgap.pdf |
| C-homog-fe | homogenized FE solve (matrix + Nitsche jump) reproduces [U]=hB to **1.6%** (−0.310 vs −0.315) | `python homog_fe.py` | REAL FE solve of the homogenized model (paper's MsFEM-Nitsche method), not a drawn step |
| FIG-panel6 | 6-window panel, BOTH panels real FE solves; jump B=−0.315 visible; B vs contrast −0.10→−0.39; oblique 2D field | `python make_panel6.py` | figures/fig_panel6.pdf (USED §5). Earlier "uniform colours" was a plotting bug (uncentred Neumann corrector + max-based colour limits); fixed by centring. |
| FIG-cellcorr | strip-cell correctors: Q⁽¹⁾ plateaus→B=−0.31; Q⁽²⁾ dipole→C/μ_m=−0.91, B₂≈0 | `python make_cellfig.py` | figures/fig_cellcorr.pdf (USED §5.2, real corrector fields) |
| C1-basis | local Neumann eigenproblem: ψ₀ constant to 1e-12 (λ₀≈0 ⇒ P¹⊂V_H); resonant modes ω₁=3820, ω₂=4380 localize in soft inclusions | `python make_basisfig.py` | figures/fig_basis.pdf — **C1 (spectral basis) computed for the first time** (was only described) |
| FIG-mesh | coarse T_H + fine resolved row + enlarged-interface (e) zoom | `python make_meshfig.py` | figures/fig_mesh.pdf (USED §3.1) |

### Figure overhaul round 2 (2026-06-17, user review of figs)
- Fig 2 (Nitsche schematic) REMOVED (duplicated eq:nitsche + the theorem).
- Fig 3 (spectral-basis schematic) REPLACED by **real computed eigenmodes** (fig_basis) — also computes C1 (closes the review gap that the spectral basis was never assembled).
- NEW mesh figure (fig_mesh, 3 views) added in §3.1 — standard in the literature; its panel (c) shows the enlarged interface thickness e explicitly.
- fig3_nitsche, fig2_spectral_basis archived to figures/_unused/. Final: 9 figures (5 real-data + domain/anatomy schematics); 13 pp, 0 undefined, 0 orphan.

### Figure overhaul round 3 (2026-06-17)
- Fig 4 (anatomy/error-decomposition flowchart, fig4_error_decomp) REMOVED (redundant "boxes-and-arrows": the 3 phases ARE the sections; the decomposition is eq:error_decomp; estimates are the Theorem). Archived to _unused/.
- FINAL FIGURE SET = 8: fig1_domain (only schematic) + 7 real-data figures (mesh, basis, cellcorr, panel6, gateref, energy, bandgap). 12 pp, 0 undefined, 0 orphan.

### Final adversarial review remediation (2026-06-17)
Final review verdict was MAJOR REVISION, driven by ONE confirmed blocker (M1, confirmed twice):
the WRITTEN variational/Nitsche form penalized [U]→0 (perfect bonding), NOT the Marigo jump, and was
inconsistent with the VALIDATED solver (homog_fe, Ventcel/spring) and with the surface mass S(ω) (absent from M_H).
FIXED by aligning the manuscript with the validated Ventcel form:
- §2.4: replaced the homogeneous Nitsche penalty by the symmetric Ventcel interface form
  P_Γ(U,V)=−(μ_m/hℬ)∫_Γ[U][V]−hC∫_Γ⟨∂_{X2}U⟩⟨∂_{X2}V⟩ (imposes the Marigo jump as natural condition;
  coercive coefficient μ_m/(h|ℬ|)>0 for ℬ<0), and added the interface surface mass hS∫_Γ⟨U⟩⟨V⟩ to the mass form m.
- §3: M_H gains the interface inertia; K_H = bulk + Ventcel P_Γ. CFL corrected to the HONEST Δt~√(Hh)/c
  (stiff physical interface spring ∝1/h), removable by implicit interface update.
- §4: a priori coercivity by sign-definiteness for ℬ<0 (M3 B>0 branch removed); energy conservation by symmetry of m,a;
  duality on the full symmetric a; scoped to ℬ<0, ℬ₂≈0 (resolves M2).
- INT-01: §5 now honestly states the convergence MMS is CONTINUOUS (isolates Term II); the jump is validated
  separately by the transmission study.
- LIT: added Pham-Maurel-Marigo (JMPS 2017, resonant) for S(ω)/band gap; Efendiev-Galvis-Hou (JCP 2013) for GMsFEM
  origin; Delourme-Haddar-Joly added to the positioning table.
- fig_bandgap recoloured (was still grayscale). Compiles 13 pp, 0 undefined, 0 bibtex warning.
REGRESSION (verified): ~18 first-review blockers/majors confirmed RESOLVED; integrity clean (no fabricated number).
Remaining nits: author affiliations placeholder; minor figure cosmetics (drift number, circle aspect).

### THIRD review (relaunch) + MATH-01 remediation (Option A) — 2026-06-17
Relaunch verdict: MAJOR REVISION, but ALL prior blockers/majors (incl. perfect-bonding M1, INT-01) confirmed
RESOLVED and integrity clean (no fabricated number). ONE new confirmed blocker MATH-01: my Ventcel rewrite was
sign-incoherent — the written coefficient −μ_m/(hℬ) is coercive but its Euler-Lagrange NEGATES the Marigo jump,
while the Marigo-faithful +μ_m/(hℬ) (= solver homog_fe κ + line 414 + eq:jump_U) is NON-coercive for ℬ<0; one
cannot have both with the bare displacement spring.
FIX (Option A, user-chosen — solver-faithful, no re-validation):
- eq:nitsche sign corrected to +μ_m/(hℬ) (= κ, matches homog_fe/line414/eq:jump_U; EL reproduces Marigo).
- Coercivity claim REPLACED: for ℬ<0 the interface term is a NEGATIVE surface stiffness ⇒ a is NOT coercive;
  well-posedness via a GÅRDING inequality (compact trace perturbation); time-domain stability via the coercive
  MODIFIED energy (Lemma a priori), not via a. "Coercive without penalty" claim removed.
- Energy theorem reframed: 𝓔_H conserved by symmetry (constant of motion) but sign-INDEFINITE for ℬ<0 (band-gap
  mechanism); boundedness from the a-priori/Gårding estimate, not from 𝓔_H positivity. Written-out energy sign fixed.
- MATH-03: Aubin–Nitsche restated as a non-conforming STRANG duality (adjoint consistency, not Galerkin orthogonality).
- MATH-04: interface consistency is O(η) (continuous P1 ⇒ residual O(h)), not O(η²); final O(η+H)/O(η+H²) UNCHANGED.
- Terminology: "symmetric Nitsche" → "symmetric Ventcel" (5×), "MsFEM--Nitsche" → "MsFEM--Ventcel", §2.4 title.
- Stale "B₂=−C₁ to 2.3%" → "to machine precision" (matches the flux-BC-corrected solve).
Compiles 13 pp, 0 undefined. Solver unchanged. Remaining: author placeholders + figure-cosmetic nits.

### FOURTH pass — relaunched review (MAJOR REVISION) + authoritative literature research — 2026-06-18
Two parallel checks: (a) relaunched 6-dimension adversarial review confirmed MATH-01 RESOLVED in the body but
re-flagged, in EVERY dimension, a surviving BLOCKER (abstract + conclusion C2 still said "coercive for the stiff
inclusion"); (b) a deep open-access literature search (PMM 2021 PMC7897651 Lemma 4.1 read verbatim;
Lombard-Maurel-Marigo JCP 2017; Darche et al. 2025; Esterhazy-Melenk; Arendt-Warma) established the AUTHORITATIVE
treatment and showed my "modified energy = E_kin+E_pot−∫[U]⟨Σ₁⟩" was mislabeled and the CFL √(Hh)/c was unsourced.
DISPUTED NUMBERS resolved by RE-RUNNING solvers (no guessing): tilted ellipse B₂=−0.04676, C₁=+0.04676 →
**machine precision CONFIRMED** (manuscript +0.047 correct; COMP-01 was reading the stale ledger row, now fixed);
band-gap edge 1.427e5/1.134e5 = **1.2584 → 1.26 CONFIRMED** (COMP-02 refuted); Bessel **5.55e-4** (text 6e-4 → 5.6e-4);
energy drift fresh run **5.33e-15** (caption → 5.3e-15). Laminate B=−0.476 vs analytic −0.477 = 0.30% (text "three
digits" → "better than 0.5%"). Circle C=−0.91407.
EDITS (text/figures only; SOLVER AND ALL PHYSICS NUMBERS UNCHANGED):
- BLOCKER: removed "and is coercive for the stiff inclusion" from abstract + conclusion C2 → honest Gårding/
  modified-energy phrasing.
- ENERGY re-anchored on the PMM 2021 enlarged-interface energy (surface kinetic carried by S + complementary
  tangential potential, C<0; the jump enters the transmission conditions, not the stored energy); the §2.4 E=½m+½a
  and the Lemma "modified energy" relabeled honestly (E sign-indefinite for B<0; the modified energy is a coercive
  auxiliary functional / Legendre transform, the analysis counterpart of the PMM positive energy). Cites PMM 2021.
- CFL reframed as the rigorous leapfrog bound Δt≤2/√(λmax(M⁻¹K)), with √(Hh)/c the resulting scaling (surface mass
  hS relaxes it, the negative stiffness tightens it); cites ChungEfendievGibson2016.
- Fully-discrete remark: "positive quadratic form" → conserved-but-indefinite E_H, with a discrete modified energy
  (adds μ_m/(h|B|)‖[U_H]‖²) supplying discrete boundedness (closes MATH-01b + COMP-04 continuous→discrete gap).
- Gårding well-posedness qualified "away from interface resonances"; cites Esterhazy-Melenk, Arendt-Warma.
- Intro "Nitsche's method" paragraph → "Weak (Ventcel) interface coupling"; dropped the false "positive-definite
  system matrix" claim (LIT-03/CF-02/W3).
- B₂=−C₁ re-attributed to Cornaggia-Touboul-Bellis 2022 (Remark 1), consistent with Marigo energy structure (LIT-02).
- S normalization reconciled: Remark 2.6 now states S=⟨ρ/ρ_m−1⟩_D=(ρ_i/ρ_m−1)|D|=0.53 (= the Sec 5.3 / Table value);
  dropped blanket "dimensionless"; added S=0.53 to Table caption (INT-S-DEFINITION/COMP-06).
- Figures regenerated with "Ventcel" titles (fig_panel6 b, fig3_energy); homog_fe.py docstring sign corrected to
  +mu_m/(h*B) (<0 for B<0); orphan fig2_convergence.pdf moved to figures/_unused/.
- references.bib: added PhamMaurelMarigo2021 (10.1098/rspa.2020.0519), CornaggiaToubouBellis2022 (10.5802/crmeca.119),
  EsterhazyMelenk2011 (10.1007/978-3-642-22061-6_9), ArendtWarma2003 (10.1023/A:1024181608863).
- Ledger: stale pre-flux-BC rows C-iface-1D/-circle/-B2C1 annotated as superseded; closure source Milton→PMM2017.
RESIDUAL HONESTY (not over-claimed in the text): the SCALAR (2017 J. Elasticity) interface-energy closed form was
NOT read verbatim (paywall); the positive-energy statement follows the PMM 2021 in-plane STRUCTURE and is hedged
accordingly. Author affiliations/emails remain placeholders (author to fill before submission).

### FIFTH pass — confirmation review (MAJOR REVISION, 0 blockers) + Lemma 4.4 proof repair — 2026-06-18
Relaunched confirmation review (54 agents): 0 blockers (the abstract/conclusion "coercive" blocker CONFIRMED gone;
all 4th-pass fixes confirmed correct AND correctly hedged; integrity clean; 2D genuine). ONE real scientific
cluster MATH-E1/E2/E3 (3 majors) — a genuine algebra error in the Lemma 4.4 (a priori) PROOF, independently
re-derived and confirmed: (E1) the interface flux in the energy identity was mis-derived ([U]<d_t Sigma_1> instead
of <Sigma_1>d_t[U]); (E2) E=1/2 m+1/2 a and the old "modified energy" could not both be conserved; (E3) the old
modified energy was built entirely from the jump term that S2.4 says is NOT stored energy, and omitted the surface
kinetic + tangential terms. FIX (user-approved, verified 3 ways: review verifier + my re-derivation + paper-internal
consistency with eq:energy_continuous + PMM structure): rewrote the Lemma proof to (i) get conservation of
E=1/2 m+1/2 a DIRECTLY from symmetry of m,a (no IBP flux detour) -- E written out with all 5 terms, sign-indefinite
for B<0; (ii) define the COERCIVE MAJORANT E_mod=E+(mu_m/h|B|)||[U]||^2>=0 (all 5 terms non-negative), explicitly
NOT the physical energy and NOT structurally the PMM energy; (iii) close the a priori bound via the compact trace
(Ehrling) embedding + Gronwall. Also: E4 "Legendre transform" -> "coercive majorant"; E5 CFL qualified to the
positive spectrum + cross-ref the discrete modified energy; NEW-1 C/mu_m units in Remark 2.6; NEW-2 typo
"B2=C1~0"->"B2,C1~0"; NEW-3 dropped B2 (out of scope) from the a-priori constant; W2 drift 5e-15->5.3e-15.
Compiles 14 pp, 0 undefined, no broken refs (removed eq:interface_perfect_derivative/modified_energy_lemma/
modified_energy_evolution labels, verified unreferenced). Solver and all physics numbers UNCHANGED; the a priori
CONCLUSION (eq:final_apriori) is unchanged. Remaining: author placeholders + data-availability URL (administrative)
and pure cosmetic nits.

### Publication styling — COLOUR + B&W-print-safe curves (2026-06-17)
Shared style module `msfem2d_verified/pubstyle.py` (SciencePlots 'science'+'no-latex' base, CM-serif mathtext,
600 dpi). Figures are FULL COLOUR (publication quality); the requirement was only that **curves remain
distinguishable when printed in black and white** — so each curve carries a distinct COLOUR **and** linestyle
**and** marker (Okabe-Ito colourblind-safe palette, pubstyle.BW). Field maps use the diverging colour map
**RdBu_r** with a thin black zero-contour (the contour also aids B&W). TikZ domain: light-blue matrix, red Γ,
blue solid jump arrow vs orange dashed mean arrow (distinguished by colour AND linestyle).
(NOTE: an earlier pass mistakenly forced everything to grayscale; reverted — colour restored, B&W-safety kept
via markers+linestyles only.) Scripts: make_meshfig/basisfig/cellfig/panel6/figs_modelA/fig_energy.py + fig1_domain.tex.

### Figure overhaul (2026-06-17) — Model A consistency + quality
- **Fig 1 (domain), Fig 3 (Nitsche), Fig 5 (anatomy)**: TikZ sources redrawn to Model A (single row→mean-line Γ, bulk=matrix, two-trace operators, γ single 1/H, strip correctors Q⁽ʲ⁾→B/B₂/S/C/C₁; no μ*, no a=e, no slab-⟨⟩_a). Each compiled + visually verified.
- **Fig 2 (scale hierarchy)**: REMOVED (redundant with nomenclature table; carried wrong "a=e").
- **Fig 4 (spectral basis)**: kept (already Model-A-consistent: enrich-not-project, ψ₀=1, P¹⊂V_H).
- **NEW Fig 5 (cellcorr)**: real strip-cell corrector fields (make_cellfig.py).
- Orphans archived to figures/_unused/. Final: 9 figures, all referenced/Model-A-consistent; 13 pp, 0 undefined, 0 orphan.

> **MANUSCRIPT INTEGRATION (2026-06-17):** main.tex fully rewritten to Model A and compiles (12 pp, 0 undefined
> refs/cites, 0 orphan figures, 0 fabricated/old-model tokens). §2 geometry+jumps+real params; §3 penalty 1/H +
> CFL H/c; §4 proofs repaired (sign-definite coercivity, full-form duality, ½a_H two cross terms, energy
> conservation split from positivity); §5 = interface params table + GATE_REF + discretization/energy + band gap;
> abstract/intro/conclusion/positioning aligned; Elsevier declarations + data-availability added; scikit-fem cited.
> Remaining: author affiliations/emails are placeholders (author to fill); optional re-run of the adversarial review.
