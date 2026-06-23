# _deprecated/

Scripts remplacés ou invalidés. Conservés pour historique uniquement.
**Ne pas utiliser pour reproduire les résultats de l'article.**

| Fichier | Raison de la dépréciation |
|---|---|
| `convergence_nitsche_mms.py` | Contient un bug d'assemblage DG (NaN sur scikit-fem < 12.0.1) ; remplacé par `src/solvers/broken_nitsche_2d.py` qui exige scikit-fem 12.0.1 |

> Voir `PROVENANCE.md` dans `data/` pour le détail des corrections apportées.
