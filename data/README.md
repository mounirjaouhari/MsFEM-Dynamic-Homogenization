# data/

Données numériques de référence et registre de provenance.

| Fichier | Contenu |
|---|---|
| `PROVENANCE.md` | Registre faisant autorité : chaque résultat numérique de l'article → commande → artéfact |
| `convergence_real.csv` | Table de convergence H-raffinement (Table 1) générée par `src/verification/regen_section5.py` |

## PROVENANCE.md

Fichier central de traçabilité. Pour chaque valeur numérique clé de l'article :
- Quelle commande la produit
- Quel script l'écrit
- Quelle tolérance est attendue

**Ce fichier NE DOIT PAS être édité manuellement.**
Il est mis à jour uniquement via `src/verification/regen_section5.py`.

## convergence_real.csv

Format :
```
H, err_L2, err_H1, rate_L2, rate_H1
```

Généré automatiquement par :
```bash
python src/verification/regen_section5.py
```

ou via :
```bash
python run_all.py --gates
```
