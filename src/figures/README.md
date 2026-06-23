# src/figures/

Scripts de génération des figures de l'article.
Chaque script produit ses figures dans `../../figures/` (dossier figures de l'article).
Tous les calculs sont réels — aucune donnée codée en dur.

| Fichier | Figures produites | Solveurs utilisés |
|---|---|---|
| `make_meshfig.py` | `fig_mesh.pdf/.png` — hiérarchie de maillages | — |
| `make_basisfig.py` | `fig_basis.pdf/.png` — fonctions de base MsFEM spectrales | — (problème local interne) |
| `make_cellfig.py` | `fig_cellcorr.pdf/.png` — correcteurs de cellule bande | `interface_params_2d` |
| `make_figs_modelA.py` | `fig_gateref`, `fig_bandgap`, `fig_field2d` (.pdf/.png) | `transmission_resolved/compare/oblique`, `resonant_mass_2d` |
| `make_panel6.py` | `fig_panel6.pdf/.png` — panneau 6 vues | `interface_params_2d`, `homog_fe`, `transmission_oblique` |
| `make_fig_energy.py` | `fig3_energy.pdf/.png` — conservation de l'énergie | `transient_energy_2d` |
| `sensitivity_wave.py` | `fig4_sensitivity.pdf/.png`, `fig5_wave_profiles.pdf/.png` | `cell_problem_2d` |

## Chemin de sortie

```python
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures")
#       src/figures/ → src/ → msfem2d_verified/ → Article 19/ → figures/
```

## Exécution

```bash
# Depuis la racine du dépôt :
python src/figures/make_cellfig.py
python src/figures/make_figs_modelA.py
# ... ou via run_all.py pour l'ordre correct
python run_all.py
```

> **Note** : exécuter via `run_all.py` garantit que les gates sont passés
> (GATE_VNV) avant la génération des figures.
