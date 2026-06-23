# src/solvers/

Solveurs par éléments finis constituant le coeur du calcul MsFEM.
Tous les scripts sont autonomes (`if __name__ == "__main__"`) et peuvent
s'exécuter directement pour vérification.

| Fichier | Rôle | Dépendances locales |
|---|---|---|
| `broken_nitsche_2d.py` | Schéma Nitsche brisé-Γ ; MMS O(H²)/O(H) ; base de toutes les tables | — |
| `cell_problem_2d.py` | Problème de cellule périodique ; μ* dans la bande Hashin-Shtrikman | — |
| `interface_params_2d.py` | Paramètres d'interface B, B₂, C, C₁, S (Table 2) via cellule bande | — |
| `homog_fe.py` | EF homogénéisé statique ; projection sur V_H Nitsche | `interface_params_2d` |
| `resonant_mass_2d.py` | Masse résonante ; carte bande interdite (Fig. 9) | — |
| `transient_energy_2d.py` | Solveur transitoire leapfrog ; dérive énergie < 1e-14 | — |
| `transmission_resolved.py` | Référence résolue : transmission onde à travers la rangée | `interface_params_2d` |
| `transmission_compare.py` | Comparaison résolu vs homogénéisé | `transmission_resolved` |
| `transmission_oblique.py` | Incidence oblique 45° | `transmission_resolved`, `interface_params_2d` |

## Exécution directe

```bash
# Vérification rapide du solveur Nitsche
python src/solvers/broken_nitsche_2d.py

# Paramètres d'interface (Table 2)
python src/solvers/interface_params_2d.py
```

## Constantes physiques partagées

| Constante | Valeur | Description |
|---|---|---|
| `MU_M` | 1.0 | Module de cisaillement matrice |
| `MU_I` | 10.0 | Module de cisaillement inclusion (régime contraste) |
| `RHO_M` | 1.0 | Densité matrice |
| `RHO_I` | 1.0 | Densité inclusion |
| `R` | 0.2821 | Rayon cercle (fraction volumique 0.25) |
| `VF` | 0.25 | Fraction volumique |
