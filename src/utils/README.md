# src/utils/

Utilitaires partagés importés par tous les scripts de génération de figures.

## pubstyle.py

Module de style matplotlib garantissant la cohérence visuelle de toutes les
figures de l'article (polices, couleurs, épaisseurs de trait, taille de page).

**Usage :**
```python
import pubstyle
pubstyle.apply()       # applique le rcParams global
ax.plot(x, y, **pubstyle.BW[0])   # style noir/blanc cycle 0
cmap = pubstyle.FIELD_CMAP        # colormap champs scalaires
```

**Ce que `apply()` configure :**
- Police : mathtext LaTeX, taille 8 pt (compatible colonnes Springer)
- Taille figure : 3.4 in (colonne simple) ou 7.0 in (pleine page)
- Lignes : lw=1.0, cycle noir/blanc avec marqueurs distincts
- Axes : pas de boîte supérieure/droite (style publication)

**Pas d'importation tierce** : uniquement `matplotlib`.
