# src/verification/

Gates de vérification numérique et régénération des artéfacts Section 5.
Ces scripts constituent la **première ligne de confiance** du dépôt :
ils doivent tous PASSER avant toute publication de figures.

| Fichier | Gate | Ce qui est vérifié |
|---|---|---|
| `smoke_poisson_mms.py` | `GATE_SMOKE` | Assemblage EF de base ; erreur L² Poisson en O(H²) |
| `nitsche_coercivity_2d.py` | `GATE_VNV[Nitsche]` | Symétrie et coercivité du schéma Nitsche brisé |
| `regen_section5.py` | `GATE_COMPUTE` | Régénère Table 1 + Table 4 ; écrit `data/convergence_real.csv` et `fig2_convergence.pdf` |

> `smoke_poisson_mms.py` et `nitsche_coercivity_2d.py` sont sans effets de
> bord (pas d'écriture fichier) — idéaux pour CI.

## Exécution

```bash
# Tous les gates seulement (rapide, ~2 min)
python run_all.py --gates

# Un gate individuel
python src/verification/smoke_poisson_mms.py
python src/verification/nitsche_coercivity_2d.py
python src/verification/regen_section5.py
```

## Sorties de regen_section5.py

| Artéfact | Destination |
|---|---|
| `convergence_real.csv` | `data/convergence_real.csv` |
| `fig2_convergence.pdf/.png` | `../figures/` (dossier article) |
| Tableau LaTeX Table 4 | `stdout` |

## Imports

`regen_section5.py` importe `broken_nitsche_2d` depuis `src/solvers/`
via le bloc `sys.path` standard en tête de fichier.
