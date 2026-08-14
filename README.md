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

### Primitive vs. conventional cell (mission #2)

`analysis/REPORT_conventional_pilot.md`: tested the pilot report's
hypothesis that the primitive/conventional cell choice explains the
near-zero correlation above, on the 6 pilot compounds. Verdict: the choice
does change the percolation weight substantially for compounds whose cell
actually differs (3x-55x, well past DFT noise), but the strongest bond
still never participates in the winning path, and the weight gets
*smaller* rather than larger as hypothesized — a bigger cell just gives
the minimum-weight search more long-range weak bonds to exploit. Extension
to the full dataset was not pursued based on this result.

### Network dimensionality + periodic min-cut (mission #3)

Two new descriptors, `percolation_path.py` untouched:
`network_dimensionality.py` (0D-3D classification via a relative bond-
strength threshold + BFS connectivity, since the existing graph's ~6 Å
cutoff makes almost everything look 3D) and `periodic_mincut.py` (minimum
total bond strength separating the crystal into two halves along a
direction, via `networkx` max-flow/min-cut on a finite ribbon graph — a
different physical question from the percolation weight: separability, not
traversability). Both validated on synthetic cases with known answers
(`tests/test_network_dimensionality.py`, `tests/test_periodic_mincut.py`)
before running on real data. See
**[`analysis/REPORT_dimensionality_mincut.md`](analysis/REPORT_dimensionality_mincut.md)**:
tested against `formation_energy_per_atom` (`mp_dataset/
fetch_formation_energy.py`) rather than `energy_above_hull`. Headline:
min-cut (normalized) is the first descriptor in this project to reach a
significant *global* correlation (ρ=0.285, p=0.0001, n=186) — stronger
than the original percolation weight ever achieved — though likely driven
partly by between-bond-type clustering rather than a clean within-type
relationship; dimensionality alone does not separate formation energy
(Kruskal-Wallis p=0.19).

### Antibonding population near the frontier (E_F/VBM)

A third, distinct question from everything above: not the integrated
ICOHP/ICOBI (a single number per bond), but *how COHP is distributed in
energy* — specifically, whether the highest-energy occupied states carry
antibonding character, by analogy with Peierls/Jahn-Teller electronic
instabilities. New module `cohp_extraction.py` (`percolation_path.py`
still untouched), built on `pymatgen.io.lobster.outputs.Cohpcar` /
`pymatgen.electronic_structure.cohp.CompleteCohp` (no hand-rolled COHPCAR
parsing). Cross-validated against the already-validated `ICOHPLIST.lobster`
across 558 bond labels in the 6 pilots (exact match for 5/6, 1e-5 eV for
the 6th); metal/gap classification cross-checked against Materials
Project rather than derived locally, which caught a real pitfall (our
LOBSTER-oriented coarse k-mesh spuriously suggests small gaps for two
known metals). See
**[`analysis/REPORT_cohp_feasibility.md`](analysis/REPORT_cohp_feasibility.md)**
(extraction pipeline validation) and
**[`analysis/METRIC_DEFINITION_antibonding.md`](analysis/METRIC_DEFINITION_antibonding.md)**
(the window/metric definition itself: one-sided window below E_F/VBM,
integrated antibonding-only COHP, raw + normalized). Validated on the 6
pilots only (synthetic numerical tests + real-data sanity checks,
`tests/test_cohp_extraction.py`) — extension to the full 186-compound
dataset, and testing whether the metric predicts anything at all, are
separate, not-yet-authorized next steps.

## Limites connues

- Suppose que `ICOHPLIST.lobster`/`ICOBILIST.lobster` proviennent du même
  calcul LOBSTER que la structure fournie (même ordre d'atomes, même
  maille) ; aucune vérification croisée de cohérence structure/POSCAR
  LOBSTER au-delà de la validation des indices d'atome.
- La complexité du Dijkstra par direction est `O(atomes)` exécutions sur un
  espace d'état de taille `atomes × (2·coord_bound+1)³` ; adapté à des
  mailles de quelques dizaines à quelques centaines d'atomes par composé,
  pas à des mailles géantes.
