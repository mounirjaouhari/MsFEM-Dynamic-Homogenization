# figures/

Figures générées automatiquement par les scripts de `src/figures/` et `src/verification/`.
Toutes les figures sont reproductibles depuis le code source via `python run_all.py`.

| Fichier | Script source | Figure article |
|---|---|---|
| `fig_mesh.pdf/.png` | `src/figures/make_meshfig.py` | Fig. 2 — hiérarchie de maillages |
| `fig_basis.pdf/.png` | `src/figures/make_basisfig.py` | Fig. 3 — fonctions de base MsFEM |
| `fig_cellcorr.pdf/.png` | `src/figures/make_cellfig.py` | Fig. 4 — correcteurs de cellule bande |
| `fig_gateref.pdf/.png` | `src/figures/make_figs_modelA.py` | Fig. 5a — validation transmission |
| `fig_bandgap.pdf/.png` | `src/figures/make_figs_modelA.py` | Fig. 5b — bande interdite |
| `fig_field2d.pdf/.png` | `src/figures/make_figs_modelA.py` | Fig. 5c — champ 2D |
| `fig_panel6.pdf/.png` | `src/figures/make_panel6.py` | Fig. 6 — panneau 6 vues |
| `fig2_convergence.pdf` | `src/verification/regen_section5.py` | Fig. convergence Table 1 |
| `fig3_energy.pdf/.png` | `src/figures/make_fig_energy.py` | Fig. 8 — conservation énergie |
| `fig4_sensitivity.pdf` | `src/figures/sensitivity_wave.py` | sensibilité μ* |
| `fig5_wave_profiles.pdf` | `src/figures/sensitivity_wave.py` | profils d'onde |

## Régénérer toutes les figures

```bash
python run_all.py
```

Les figures PDF sont incluses dans le dépôt pour permettre une inspection directe
sans avoir à exécuter le code. Les calculs (≥ 10 min) ne sont nécessaires que pour
vérifier la reproductibilité ou modifier les paramètres.
