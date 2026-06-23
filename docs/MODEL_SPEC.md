# MODEL SPEC — Model A: Single Row of Elastic Inclusions

**Version:** 2026-06-16 (locked, supersedes all earlier model descriptions)
**Referenced by:** `data/PROVENANCE.md` (5 occurrences), all solver scripts

This document defines the mathematical model implemented in this repository.
It is the authoritative source for all physical constants, corrector definitions,
and parameter formulas used in the paper.

---

## 1. Physical configuration

A single periodic row of elastic inclusions is embedded in a homogeneous
matrix. The row is confined to a thin layer of thickness `e = O(h)` about
the mean line `X₁ = 0`.

**Bulk subdomains:**
- `Ω⁻ = {X₁ < -e/2}` — pure matrix (μ_m, ρ_m)
- `Ω⁺ = {X₁ > +e/2}` — pure matrix (μ_m, ρ_m)

**No bulk homogenization:** A single row does not produce a bulk effective
medium. All microstructural effects are captured by the effective interface
jump conditions across `Γ = {X₁ = 0}`.

---

## 2. Material parameters

### Structural regime (stiff inclusion)
| Parameter | Symbol | Value |
|-----------|--------|-------|
| Matrix shear modulus | μ_m | 12.0 GPa |
| Matrix density | ρ_m | 2500 kg/m³ |
| Inclusion shear modulus | μ_i | 78.0 GPa (contrast μ_i/μ_m = 6.5) |
| Inclusion density | ρ_i | 7800 kg/m³ (contrast ρ_i/ρ_m = 3.1) |
| Inclusion radius | R | 0.282 (in units of period h = 1) |
| Area fraction | vf | πR² ≈ 0.25 |

### Locally-resonant regime (soft inclusion)
| Parameter | Symbol | Value |
|-----------|--------|-------|
| Matrix shear modulus | μ_m | 12.0 GPa |
| Matrix density | ρ_m | 2500 kg/m³ |
| Inclusion shear modulus | μ_i | 0.5 GPa (contrast μ_i/μ_m = 0.042) |
| Inclusion density | ρ_i | 9000 kg/m³ |
| Inclusion radius | R | 0.282 |

---

## 3. Governing equation

Anti-plane shear wave equation:

```
ρ(X) ∂ₜₜU − ∇·(μ(X) ∇U) = f(X, t)    in Ω \ Γ
```

with `ρ(X)`, `μ(X)` piecewise constant: matrix values in `Ω±`, inclusion
values inside the layer `|X₁| ≤ e/2`.

---

## 4. Effective interface (homogenized model)

After matched-asymptotic analysis (Marigo–Maurel–Pham–Sbitti, J. Elasticity 2017),
the inclusion row is replaced by jump conditions at `Γ = {X₁ = 0}`:

```
[U]     = h ( B ⟨∂_{X₁}U⟩ + B₂ ∂_{X₂}⟨U⟩ )
[Σ₁]    = h ( C ∂²_{X₂}⟨U⟩ + ρ_m S ω² ⟨U⟩ )
```

where `[·]` is the jump across `Γ`, `⟨·⟩` is the symmetric mean,
`Σ₁ = μ_m ∂_{X₁}U`, and `S` is the effective surface mass.

The four dimensionless interface parameters `B, B₂, C, C₁` and the static
surface mass `S` are defined and computed in §6 below.

---

## 5. Weak (Ventcel) formulation

The interface conditions are imposed weakly. The bilinear forms are:

**Mass form:**
```
m(U, V) = ∫_Ω ρ_m UV + h S ∫_Γ ⟨U⟩⟨V⟩
```

**Stiffness form:**
```
a(U, V) = ∫_Ω μ_m ∇U·∇V + P_Γ(U, V)
```

where the Ventcel interface form is:

```
P_Γ(U, V) = (μ_m/h) ∫_Γ [ -B [U]⟨∂_{X₁}V⟩ - B [V]⟨∂_{X₁}U⟩
                            + (μ_m/|B|) [U][V] ]
           + h C ∫_Γ ∂_{X₂}⟨U⟩ ∂_{X₂}⟨V⟩
```

The penalty coefficient `μ_m / (h |B|)` is the physical surface spring
(scales as O(1/h)); it sets the largest eigenvalue `λ_max ~ c²/(Hh)` and
the CFL condition `Δt ≤ C_CFL √(Hh) / c`.

---

## 6. Strip-cell corrector problems and interface parameters

### Geometry
Infinite strip `Y_∞ = ℝ_{y₁} × (0,1)`, 1-periodic in `y₂`, truncated to
`(-L, L) × (0,1)` with `L = 10` (far-field plateau reached).
One circular inclusion of radius `R` centred at `(0, 1/2)`.

### Normal corrector Q¹
Solve:
```
∇·(μ(y) ∇u) = 0    in (-L, L) × (0,1)
u = ±L              at y₁ = ±L  (matrix far field)
1-periodic in y₂
```
Define `Q¹ = u − y₁`.

### Tangential corrector Q²
Solve:
```
∇·(μ(y)(∇Q² + e₂)) = 0    in (-L, L) × (0,1)
∂_{y₁} Q² = 0              at y₁ = ±L  (Neumann)
1-periodic in y₂
```

### Interface parameters (dimensionless, multiply by period h in jump law)
```
B  = ⟨Q¹⟩_{y₁=+x*} − ⟨Q¹⟩_{y₁=−x*}     (normal displacement-jump compliance)
B₂ = ⟨Q²⟩_{y₁=+x*} − ⟨Q²⟩_{y₁=−x*}     (tangential compliance)
C  = ∫ μ ∂_{y₂} Q² dy                     (tangential stress coupling)
C₁ = ∫ μ ∂_{y₂} Q¹ dy                     (mixed coupling)
S  = (ρ_i − ρ_m) · π R²                   (static excess surface mass)
```

where `x*` is taken beyond the inclusion (plateau region, `x* = 3R`).

### Energy-consistency identity (Cornaggia–Touboul–Bellis, CRAS 2022, Rem. 1)
```
B₂ = −C₁
```
This relation is verified to machine precision for all configurations in
`interface_params_2d.py`.

### Verified values (structural regime)

| Configuration | B | B₂ | C/μ_m | C₁/μ_m | S |
|---|---|---|---|---|---|
| Matrix (control) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 1D laminate | −0.476 | 0.000 | 0.000 | 0.000 | 0.000 |
| Centred circle | −0.311 | ≈0 | −0.914 | ≈0 | 0.530 |
| Tilted ellipse | −0.18 | −0.047 | −0.44 | +0.047 | — |

Source: `python interface_params_2d.py` (logged in `data/PROVENANCE.md`).

---

## 7. Resonant surface mass (locally-resonant regime)

The static `S` is replaced by a frequency-dependent surface mass:

```
S(ω) = S₀ + Σ_n  aₙ ωₙ² / (ωₙ² − ω²)
```

where `ωₙ` are the eigenfrequencies of the clamped inclusion and
`aₙ = 2ωₙ⁻¹ ∫_Y ρ (ψₙ)² > 0` are the modal oscillator strengths
of the n-th symmetric (breathing) mode `ψₙ`.

Only symmetric (m=0) breathing modes contribute; dipole modes
(m=1) have `aₙ ≈ 0` (verified in `resonant_mass_2d.py`).

`S(ω)` changes sign above each resonance `ωₙ`, opening a
**negative-effective-mass band gap** (metamaterial mechanism).

Verified: FE ω₁ matches analytic clamped-disk value `j₀₁ cᵢ/R` to `5.6e-4`.

---

## 8. What this model is NOT

- **Not** a bulk periodic composite homogenization (μ* ≠ μ_m in bulk).
  The old `cell_problem_2d.py` result (μ* = 17.47 GPa) is used only as the
  Nitsche penalty scale — it is NOT the bulk modulus of the homogenized medium.
- **Not** a two-scale bulk effective medium. The Bensoussan–Lions–Papanicolaou
  theory does not apply to a single-row geometry.
- **Not** valid for multiple rows or bulk periodic arrays.

---

## 9. Sources

- Marigo, Maurel, Pham, Sbitti — *J. Elasticity* 128 (2017) 265–289
- Pham, Maurel, Marigo — *JMPS* 106 (2017) 80–94
- Cornaggia, Touboul, Bellis — *C. R. Mécanique* 350 (2022)
- Lombard, Maurel, Marigo (preprint, arXiv)
