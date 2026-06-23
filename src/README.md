# src/

Code source Python du dépôt MsFEM-Dynamic-Homogenization.

| Sous-dossier | Rôle |
|---|---|
| [`solvers/`](solvers/) | Solveurs EF (Nitsche brisé, problèmes de cellule, énergie, transmission) |
| [`figures/`](figures/) | Scripts de génération des figures de l'article |
| [`verification/`](verification/) | Gates de vérification numérique et régénération Section 5 |
| [`utils/`](utils/) | Utilitaires partagés (style matplotlib) |

## Dépendances entre sous-dossiers

```
utils/pubstyle.py
    ↑ importé par figures/*

solvers/{broken_nitsche_2d, cell_problem_2d, interface_params_2d, ...}
    ↑ importés par figures/* et verification/*

verification/regen_section5.py
    → importe broken_nitsche_2d (solvers/)
    → écrit data/convergence_real.csv et ../figures/*.pdf

figures/make_*.py, sensitivity_wave.py
    → importent pubstyle (utils/) + solveurs (solvers/)
    → écrivent ../../figures/*.pdf (dossier figures de l'article)
```

## Résolution des imports

Chaque script dans `figures/`, `verification/` commence par un bloc :

```python
import sys as _sys, os as _os
_SRC = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for _d in ('solvers', 'utils'):
    _p = _os.path.join(_SRC, _d)
    if _p not in _sys.path: _sys.path.insert(0, _p)
del _sys, _os, _SRC, _d, _p
```

Cela rend `solvers/` et `utils/` visibles sans modifier `PYTHONPATH`.
