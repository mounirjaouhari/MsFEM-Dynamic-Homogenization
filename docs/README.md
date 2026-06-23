# docs/

Documentation mathématique et spécifications du modèle.

| Fichier | Contenu |
|---|---|
| `MODEL_SPEC.md` | Spécification complète du modèle : géométrie, paramètres matériaux, équation gouvernante, conditions d'interface, formulation Ventcel, correcteurs de cellule bande |

## MODEL_SPEC.md

Document de référence pour :
- La géométrie exacte (rangée périodique de fibres circulaires, fraction volumique 0.25)
- Les deux régimes : contraste (μ_I/μ_M = 10) et résonant (μ_I/μ_M = 100)
- Les conditions de saut d'interface exactes vérifiées
- Les valeurs numériques de B, B₂, C, C₁, S reproduites par `src/solvers/interface_params_2d.py`
- Les formules analytiques de vérification (Hashin-Shtrikman, masse de Drude-Lorentz)

Ce fichier est la source de vérité pour tout nouveau développeur qui veut
comprendre CE QUE calcule le code avant de l'exécuter.
