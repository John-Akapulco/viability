# percolation_path.py

Descripteur de « chemin de percolation de moindre résistance » calculé par
post-traitement sur une structure cristalline relaxée et ses données
ICOHP/ICOBI (LOBSTER). Aucune dynamique moléculaire, aucun calcul NEB,
aucune supercell physique : la périodicité est gérée par étiquetage
vectoriel des arêtes du graphe de liaisons (voltage graph / gain graph),
et le chemin non contractile de poids minimal est trouvé par un Dijkstra
sur l'espace d'état étendu `(atome, translation cumulée)`.

## Installation

```bash
pip install -r requirements.txt   # pymatgen
```

Nécessite Python ≥ 3.9 (utilise `pymatgen.io.lobster.outputs.Icohplist`).

## Format d'entrée attendu

```
dataset/
  compound_A/
    POSCAR (ou CONTCAR, ou *.cif)      # structure relaxée, maille primitive
    ICOHPLIST.lobster                   # sortie LOBSTER, avec colonnes de translation
    ICOBILIST.lobster                   # optionnel
  compound_B/
    ...
```

Chaque sous-répertoire de `--root` est traité comme un composé indépendant.
Le fichier `ICOHPLIST.lobster`/`ICOBILIST.lobster` doit inclure la colonne
de translation de réseau (`tx ty tz`) par liaison — c'est le cas standard
depuis LOBSTER ≥ 3. C'est cette information qui permet de reconstruire le
graphe périodique étiqueté sans jamais dupliquer physiquement les atomes.

## Utilisation

```bash
python percolation_path.py --root dataset --metric icohp --output results.csv
python percolation_path.py --root dataset --metric icohp,icobi --output results.csv --also-json results.json
python percolation_path.py --root dataset --metric icohp --bond-pair Fe-O --output results.csv
```

Options principales :
- `--metric icohp|icobi|icohp,icobi` : quelle(s) grandeur(s) utiliser comme
  poids d'arête (|ICOHP| ou |ICOBI|). Les deux peuvent être calculées côte
  à côte pour comparaison.
- `--bond-pair Fe-O` : restreint le graphe et les agrégats aux liaisons
  entre ces deux espèces (sinon toutes les liaisons du fichier sont
  utilisées).
- `--coord-bound N` : borne du domaine d'exploration des translations
  cumulées pendant la recherche de chemin. Par défaut, dérivée
  automatiquement du plus grand vecteur de translation présent dans les
  données d'entrée (donc directement liée au rayon de coupure du calcul
  ICOHP/ICOBI amont, pas un paramètre de convergence à ajuster).

## Sortie

Un enregistrement par composé (CSV, une ligne par composé ; ou JSON,
structure imbriquée avec le détail par direction). Colonnes principales
(préfixées par le nom de la métrique, ex. `icohp_...`) :

- `*_percolation_weight_min` : poids du chemin de percolation le plus
  faible, toutes directions confondues — le descripteur principal.
- `*_percolation_direction` : direction (`a`, `b`, ou `c`) correspondante.
- `*_percolation_weight_a/b/c` et `*_percolation_status_a/b/c` : poids
  minimal par direction, avec statut explicite (`ok` ou `disconnected` si
  aucun chemin non contractile n'existe dans cette direction — jamais une
  valeur infinie silencieuse).
- `*_sum`, `*_mean`, `*_min`, `*_max` : agrégats classiques sur les mêmes
  liaisons, pour comparaison directe avec le descripteur de percolation.
- `error`, `warnings` : diagnostics par composé (le traitement par lot ne
  s'interrompt jamais sur un composé en échec — l'erreur est consignée
  dans la ligne correspondante et le lot continue).

## Algorithme (résumé)

1. Chaque liaison ICOHP/ICOBI devient une arête `(atome_i, atome_j, (nx,ny,nz))`
   du graphe périodique, avec poids `|valeur|`. Les deux sens de l'arête
   sont ajoutés (translation opposée dans l'autre sens).
2. Pour chaque direction de réseau `a`, `b`, `c`, on cherche — pour chaque
   atome de départ possible — le chemin de poids minimal dans l'espace
   d'état `(atome, translation cumulée)` reliant `(atome, (0,0,0))` à
   `(atome, direction)`. C'est un Dijkstra classique (poids ≥ 0 car
   valeurs absolues), valide car l'espace d'état est fini une fois borné
   par `--coord-bound`.
3. Le minimum sur les atomes de départ donne le poids de percolation de
   cette direction ; le minimum sur les trois directions donne le
   descripteur principal.
4. Si aucun état-cible n'est atteint dans une direction (réseau de
   liaisons déconnecté selon cette direction compte tenu du rayon de
   coupure ICOHP/ICOBI), c'est signalé explicitement (`status:
   disconnected`), jamais silencieusement comme un poids infini ou nul.

## Tests

```bash
python -m unittest discover -s tests -v
```

Couvre : cas isotrope (poids identique dans les 3 directions), cas
anisotrope (direction la plus faible correctement identifiée), cas
déconnecté (signalé explicitement, pas de crash ni de valeur erronée),
cas où un chemin indirect à plusieurs sauts bat une liaison directe plus
coûteuse (validation que l'algorithme résout bien un plus-court-chemin
global et non une simple sélection de la liaison la plus faible), et un
test d'intégration bout-en-bout avec un `ICOHPLIST.lobster` synthétique
parsé via pymatgen.

## Exemples

`examples/dataset/` contient trois composés jouets illustrant les trois
cas ci-dessus (anisotrope, chemin indirect moins cher que la liaison
directe, direction déconnectée) :

```bash
python percolation_path.py --root examples/dataset --metric icohp --output /tmp/example_results.csv -v
```

## Analysis

`analysis/` contains a statistical study of the percolation descriptor
against thermodynamic stability (`energy_above_hull`) over the 186-compound
dataset in `mp_dataset/structures/` (60 experimental-stable + 60
experimental-metastable + 60 theoretical-metastable, from `mp_dataset/
select_campaign.py`, plus the 6-compound pilot). See
**[`analysis/REPORT.md`](analysis/REPORT.md)** for the full write-up:
Spearman correlations (overall and per bond_type) for the raw and
normalized percolation weight vs. the classic ICOHP aggregates, a reference
logistic regression (stable vs. metastable) with cross-validated AUC,
figures, compute-time table, and a documented limits/next-steps section.
Headline result: no significant overall correlation with `energy_above_hull`
at this sample size, and no clear predictive edge over the classic
aggregates in the one subgroup (metallic) that shows a nominally
significant signal — see the report for why, and what would need to change
before this justifies a symbolic-regression (SISSO) pass.

Pipeline: `analysis/build_dataset.py` (join `percolation_path.py` output +
MP metadata → `analysis/percolation_vs_hull.csv`) then
`analysis/stats_analysis.py` (correlations, logistic regression, figures →
`analysis/stats_summary.json` + `analysis/figures/`).

## Limites connues

- Suppose que `ICOHPLIST.lobster`/`ICOBILIST.lobster` proviennent du même
  calcul LOBSTER que la structure fournie (même ordre d'atomes, même
  maille) ; aucune vérification croisée de cohérence structure/POSCAR
  LOBSTER au-delà de la validation des indices d'atome.
- La complexité du Dijkstra par direction est `O(atomes)` exécutions sur un
  espace d'état de taille `atomes × (2·coord_bound+1)³` ; adapté à des
  mailles de quelques dizaines à quelques centaines d'atomes par composé,
  pas à des mailles géantes.
