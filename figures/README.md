# figures/

Figures générées automatiquement par les scripts de `src/figures/` et `src/verification/`.
Toutes les figures sont reproductibles depuis le code source via `python run_all.py`.

| Fichier | Script source | Figure article (main-sn.tex) |
|---|---|---|
| `fig1_domain.pdf` | `fig1_domain.tex` (TikZ, manuel) | Fig. 1 — domaine physique et homogénéisé |
| `fig_mesh.pdf/.png` | `src/figures/make_meshfig.py` | Fig. 2 — hiérarchie de maillages |
| `fig_basis.pdf/.png` | `src/figures/make_basisfig.py` | Fig. 3 — fonctions de base spectrales MsFEM |
| `fig_cellcorr.pdf/.png` | `src/figures/make_cellfig.py` | Fig. 4 — correcteurs de cellule bande |
| `fig_gateref.pdf/.png` | `src/figures/make_figs_modelA.py` | Fig. 5 — validation transmission (gate/ref) |
| `fig_panel6.pdf/.png` | `src/figures/make_panel6.py` | Fig. 6a — panneau 6 vues saut Ventcel |
| `fig_scattering_3times.pdf/.png` | `src/figures/make_fig_scattering.py` | Fig. 6b — diffusion transitoire DNS vs EI |
| `fig3_energy.pdf/.png` | `src/figures/make_fig_energy.py` | Fig. 7 — conservation énergie leapfrog |
| `fig_bandgap.pdf/.png` | `src/figures/make_figs_modelA.py` | Fig. 8 — bande interdite à masse effective négative |
| `fig_field2d.pdf/.png` | `src/figures/make_figs_modelA.py` | (hors article — champ 2D oblique, référence interne) |
| `fig4_sensitivity.pdf` | `src/figures/sensitivity_wave.py` | (hors article — sensibilité μ* vs contraste) |
| `fig5_wave_profiles.pdf` | `src/figures/sensitivity_wave.py` | (hors article — profils d'onde snapshot) |
| `fig2_convergence.pdf` | `src/verification/regen_section5.py` | (hors article — courbes de convergence Table 1) |
| `anim_wave.gif` | `src/figures/make_anim_wave.py` | Animation — propagation d'onde (GitHub uniquement) |
| `anim_bandgap.gif` | `src/figures/make_anim_bandgap.py` | Animation — bande interdite (GitHub uniquement) |
| `anim_scattering.gif` | `src/figures/make_anim_scattering.py` | Animation — diffusion transitoire DNS vs EI (GitHub uniquement) |

## Cohérence avec l'article (main-sn.tex, 8 figures)

`fig_panel6.pdf` et `fig_scattering_3times.pdf` partagent **un seul environnement figure**
(label `fig:panel6_and_scattering`) — ils comptent ensemble pour Fig. 6.

| Fig. | Label LaTeX | Fichier(s) généré(s) | Script | Paramètres clés |
|------|-------------|----------------------|--------|-----------------|
| 1 | `fig:domain` | `fig1_domain.pdf` | TikZ — | — |
| 2 | `fig:mesh` | `fig_mesh.pdf` | `make_meshfig.py` | — |
| 3 | `fig:spectral_basis` | `fig_basis.pdf` | `make_basisfig.py` | MU_I=0.5e9 (mou, localement résonant) |
| 4 | `fig:cellcorr` | `fig_cellcorr.pdf` | `make_cellfig.py` | μi/μm=6.5 → B=−0.31, C/μm=−0.91 |
| 5 | `fig:gateref` | `fig_gateref.pdf` | `make_figs_modelA.py` | B_eff→−0.30, S_eff→0.51 |
| 6a | `fig:panel6_and_scattering` | `fig_panel6.pdf` | `make_panel6.py` | saut [U]=hB à 1.6% |
| 6b | (même figure) | `fig_scattering_3times.pdf` | `make_fig_scattering.py` | DNS vs EI, μi/μm=6.5, ρi/ρm=3.1, R=0.282 |
| 7 | `fig:energy` | `fig3_energy.pdf` | `make_fig_energy.py` | 4000 pas, drift ~5.8×10⁻¹⁵ |
| 8 | `fig:bandgap` | `fig_bandgap.pdf` | `make_figs_modelA.py` | f₁≈18 kHz, gap [1, 1.26]ω₁ |

## Régénérer toutes les figures

```bash
python run_all.py
```

Les figures PDF sont incluses dans le dépôt pour permettre une inspection directe
sans avoir à exécuter le code. Les calculs (≥ 10 min au total) ne sont nécessaires
que pour vérifier la reproductibilité ou modifier les paramètres.
