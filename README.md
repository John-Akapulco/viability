# viability

Do local, LOBSTER-derived ICOHP/ICOBI descriptors of crystal chemical
bonding predict thermodynamic stability? Tested over a growing set of
VASP+LOBSTER calculations (`mp_dataset/structures/`, 349+ compounds),
five missions in, across four distinct descriptor families.

**Flagship result so far**: not the integrated per-bond ICOHP/ICOBI (a
single number per bond, the project's original angle — background
below), but *how COHP is distributed in energy* — specifically, the
fraction of antibonding character occupied just below the Fermi
level/VBM (`cohp_extraction.py`, mission #4). At the current 349-compound
scale it reaches ρ=−0.282, p=9.2×10⁻⁸ (n=346) against
`formation_energy_per_atom`, and its `bond_type=covalent` subgroup
(n=33, ρ=−0.551, p=0.0009) remains the most robust single-stratum result
in the project — it survived a near-50% growth in sample size and got
*stronger*, not weaker. Numerically the strongest global correlation is
now **reaction ICOHP** (mission #5, ρ=0.502, n=275) — see below — which
also newly shows a real (if not yet mechanistically understood)
`is_metal`-stratified signal, a first for that descriptor.
See [Antibonding population near the frontier](#antibonding-population-near-the-frontier-ef-vbm-mission-4) and [Reaction ICOHP](#reaction-icohp-mission-5) below.

Every descriptor here shares one methodological root: `percolation_path.py`
(the project's first descriptor, detailed near the bottom of this file)
established the periodic-graph representation of ICOHP/ICOBI bond data —
translation-labeled edges, no physical supercell duplication — and the
statistical conventions every later mission reuses without exception:
Spearman correlation only (no assumed linearity), mandatory `bond_type`/
`is_metal` stratification before any global correlation is trusted,
n<15 subgroups always flagged, and no symbolic regression (SISSO) until
a descriptor's signal is clean and well-powered enough per stratum to
justify one — still not the case for any descriptor as of mission #5.

## Descriptors (ordered by current headline strength, not by mission number)

1. **[Reaction ICOHP](#reaction-icohp-mission-5)** — mission #5, `reaction_icohp.py`. ICOHP analog of formation energy; numerically the strongest correlation in the project (ρ=0.502), and — as of the 349-compound scale — the first descriptor whose signal survives an `is_metal` stratification rather than being fully confound-explained.
2. **[Antibonding population near the frontier (E_F/VBM)](#antibonding-population-near-the-frontier-ef-vbm-mission-4)** — mission #4, `cohp_extraction.py`. Energy-resolved COHP, not integrated ICOHP. Its `bond_type=covalent` subgroup is the most robust single-stratum result in the project.
3. **[Network dimensionality + periodic min-cut](#network-dimensionality--periodic-min-cut-mission-3)** — mission #3, `network_dimensionality.py` + `periodic_mincut.py`. Graph separability, not traversability.
4. **[percolation_path.py — the original descriptor](#percolation_pathpy--the-original-descriptor-mission-1)** — mission #1. Minimum-weight non-contractile path through the periodic bond graph. No significant correlation found on its own, but the methodological foundation (periodic graph construction, statistical conventions) every descriptor above still builds on.

---

## Antibonding population near the frontier (E_F/VBM) (mission #4)

A distinct question from integrated ICOHP/ICOBI (a single number per
bond): not *how much* bonding a compound has in total, but *how COHP is
distributed in energy* — specifically, whether the highest-energy
occupied states carry antibonding character, by analogy with
Peierls/Jahn-Teller electronic instabilities. `cohp_extraction.py`, built
on `pymatgen.io.lobster.outputs.Cohpcar` / `pymatgen.electronic_structure.cohp.CompleteCohp`
(no hand-rolled COHPCAR parsing). Cross-validated against the
already-validated `ICOHPLIST.lobster` across 558 bond labels in the 6
pilot compounds (exact match for 5/6, 1e-5 eV for the 6th); metal/gap
classification cross-checked against Materials Project rather than
derived locally, which caught a real pitfall (our LOBSTER-oriented coarse
k-mesh spuriously suggests small gaps for two known metals). See
**[`analysis/REPORT_cohp_feasibility.md`](analysis/REPORT_cohp_feasibility.md)**
(extraction pipeline validation) and
**[`analysis/METRIC_DEFINITION_antibonding.md`](analysis/METRIC_DEFINITION_antibonding.md)**
(the window/metric definition itself: one-sided window below E_F/VBM,
integrated antibonding-only COHP, raw + normalized). Validated on the 6
pilots only (synthetic numerical tests + real-data sanity checks,
`tests/test_cohp_extraction.py`) before any extension.

Extended to the full 186-compound dataset in `analysis/compute_antibonding_all.py`
and tested against `formation_energy_per_atom` in
`analysis/stats_analysis_antibonding.py`; both rerun unmodified over the
349-compound dataset once the `extension4` batch (89 alkali/alkaline-earth
binaries) was added (347/349 succeeded, 2 pre-existing COD-sourced gaps).
See **[`analysis/REPORT_antibonding.md`](analysis/REPORT_antibonding.md)**.
Headline: the normalized metric reaches ρ=−0.282, p=9.2×10⁻⁸ (n=346) —
weaker than the ρ=−0.328 (n=186) originally reported, but still comfortably
ahead of every other descriptor except reaction ICOHP (min-cut's own
headline, in fact, collapsed to non-significance at this scale — see the
report). The sign (more antibonding population near the frontier
associates with *more negative*, i.e. more stable, formation energy) is
still the opposite of the naive Peierls/Jahn-Teller reading, and
diagnostics still show the global number is substantially a between-group
effect — but `is_metal=False` now also survives on its own (n=157,
p=0.001), which it did not at n=186. **`bond_type=covalent` (n=33) not
only still holds up under stratification but got *stronger* with more
data** (ρ=−0.551, p=0.0009, vs. the original ρ=−0.490 at n=23) — still
the single most robust bond-type-stratified result in the project. The
old `bond_type=ionic` result (n=9, ρ=−0.683) did not replicate once the
sample grew to n=53 (ρ=−0.184, not significant) — a clean demonstration
of why n<15 rows are flagged rather than trusted. See the report for the
full within-group diagnostics, the ΔE-sensitivity check (still robust
across 0.5–2.0 eV), and why the sign is not yet interpretable as
confirming or refuting the instability hypothesis.

`ICOBI`-based near-frontier windowing (as opposed to `ICOHP`/COHP) is a
natural extension not yet implemented — `percolation_path.py` already
treats ICOHP and ICOBI symmetrically as alternative edge weights, but
`cohp_extraction.py` is COHP/ICOHP-only for now.

## Reaction ICOHP (mission #5)

A thermochemistry-flavored question: not a compound's own bonding
topology or energy distribution in isolation, but whether its total ICOHP
is "worth more" than the same atoms would have in a reference
configuration — the ICOHP analog of `formation_energy_per_atom`. New
module `reaction_icohp.py`, three reaction types (decomposition into
elements, polymorph comparison, decomposition into a compound +
elements), balanced via `pymatgen.analysis.reaction_calculator.Reaction`.
Defined and validated on n=1–5 hand-worked real examples (Ca3N2, Mn2O7,
carbon allotropes, TiO2 high-pressure polymorphs) in
**[`analysis/METRIC_DEFINITION_reaction_icohp.md`](analysis/METRIC_DEFINITION_reaction_icohp.md)**,
which already flagged the key caveat before any statistics were run:
ICOHP sees orbital-overlap bond population, not electrostatic/Madelung or
van der Waals energy, so strongly ionic compounds and van-der-Waals-bound
polymorphs (e.g. graphite) are expected to misbehave.

Case 1 (decomposition into elements), using 62 elemental reference
calculations (`mp_dataset/download_elements_reference.py` + hand-picked
extension compounds), extended from 192 to **281 compounds** once
`extension4` (89 alkali/alkaline-earth binaries) was added — all 89 new
compounds succeeded with 0 missing elemental references, since
`extension4`'s chemistry was already fully covered. Tested in
`analysis/stats_analysis_reaction_icohp.py`. See
**[`analysis/REPORT_reaction_icohp.md`](analysis/REPORT_reaction_icohp.md)**.
Headline: ρ=0.502, p=6.3×10⁻¹⁹ (n=275) against `formation_energy_per_atom`
— up from ρ≈0.37 at n=186, and now the clear strongest global correlation
in the project by a wider margin than before. The `bond_type`-level
between-group effect (metal vs. gapped, and across the three `bond_type`
strata) still explains that layer of the pooled signal entirely — no
`bond_type` stratum reaches significance on its own, same as before — but
**the coarser `is_metal` split now tells a different story**: both
`is_metal=False` (ρ=0.319, p=0.0001) and `is_metal=True` (ρ=0.256,
p=0.0025) newly survive, where at n=186 neither did. This is the first
descriptor in the project whose signal is not fully explained away by
between-group clustering at every level tested. Against `energy_above_hull`
the correlation is no longer flatly null either (ρ=0.14, p=0.018, still an
order of magnitude weaker than the formation-energy result). Case 2
(polymorph comparison), extended from 8 to **49 polymorph groups**
specifically because `extension4` was designed to fix this project's
polymorph-density gap: 53.1% agreement between most-bonding and
most-stable member, statistically indistinguishable from the 45.9% chance
baseline (P=0.19) — a more informative null result than the old 8-group
pass, still concluding bonding does not track polymorph stability. Case 3
(decomposition into a compound + elements) remains not attempted; see the
report for why.

**`reaction_analysis/` (new, schema-driven redesign of this axis)**:
a from-scratch Pydantic schema (`CompoundEntry`, `Reaction`,
`ReactionResult`) meant to eventually cover all three reaction types
above through one common, testable data model, rather than the ad hoc
`reaction_icohp.py` functions above. Ships with `parse_lobster.py`
(builds a `CompoundEntry` from `ICOHPLIST.lobster`/`ICOBILIST.lobster` +
structure, with an explicit regression test confirming LOBSTER lists each
periodic bond once, not once per direction — the assumption
`sum_total_eV`'s unfiltered summation depends on), `balance.py`
(element-by-element stoichiometric balance checking, coefficient
auto-derivation for decomposition-into-elements), and `delta.py` (the
three ΔICOHP/ΔICOBI normalizations — per formula unit, per atom, and a
non-conservative per-bond diagnostic — computed together, never one in
isolation). Populated with a first real-data batch (6 compounds, 3
`decomposition_to_elements` reactions), cross-checked 3/3 against
`reaction_icohp.py`'s numbers, then extended to the full case-1 population
(`analysis/populate_reaction_analysis_case1_full.py`, `reactions_dataset/`)
— 192/192, then 281/281 once `extension4` was added, match
`reaction_icohp_case1.csv` after the known sign flip. Its sign convention
(products − reactants) is the *opposite* of `reaction_icohp.py`'s — the
two are not interchangeable.

**Endobondic/exobondic classification** (Reitz & Dronskowski,
ic-2026-04181q): `nearest_neighbor.py` (first-coordination-shell bond
filtering — a relative, self-calibrating gap detector, no hardcoded
Angstrom cutoff), `classify.py` (`BondingLabel`
endobondic/exobondic from ΔICOHP's sign, and a deliberately more cautious
`ViabilityLabel` — `UNSTABLE_NONEXISTENT` always carries a warning that an
exobondic sign never proves non-existence, per the manuscript's own
Mn2O7 caveat), and `units.py` (eV ↔ kJ/mol). `parse_lobster.py` gained an
opt-in `bond_filter="nearest_neighbor"` (default stays `"unfiltered"`, to
avoid silently changing the already-validated case-1 population above),
and `delta.py`'s `ReactionResult` now carries a `bonding_label`. Validated
against all 7 worked examples from the manuscript
(`tests/test_reitz_dronskowski_validation.py`,
`tests/fixtures/reitz_dronskowski_cases.yaml`) — published, external
numbers, independent of this project's own CSP data. **Not yet connected
to real LOBSTER data for classification** (the 7-case validation uses
literature bond values directly, bypassing `parse_lobster.py`) — that
remains a follow-up step.

## Network dimensionality + periodic min-cut (mission #3)

Two descriptors, `percolation_path.py` untouched:
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
min-cut (normalized) was the first descriptor in the project to reach a
significant *global* correlation (ρ=0.285, p=0.0001, n=186) — stronger
than the original percolation weight ever achieved — though likely driven
partly by between-bond-type clustering rather than a clean within-type
relationship; dimensionality alone does not separate formation energy
(Kruskal-Wallis p=0.19).

---

## `percolation_path.py` — the original descriptor (mission #1)

The project's first descriptor and the reason everything above shares a
common data model. Post-processing on a relaxed crystal structure and its
ICOHP/ICOBI (LOBSTER) data: no molecular dynamics, no NEB calculation, no
physical supercell — periodicity is handled by vector-labeling the bond
graph's edges (a voltage graph / gain graph), and the minimum-weight
non-contractile path is found by Dijkstra search over the extended state
space `(atom, cumulative translation)`.

Tested against `energy_above_hull` (mission #1) and, on the primitive-vs-
conventional-cell question, revisited on the 6 pilot compounds (mission
#2, see below): **no significant correlation was found on its own** — see
[Analysis](#analysis-of-percolation_pathpy-missions-1-2) — but the graph
construction and Dijkstra machinery here are exactly what every later
descriptor's own graph (min-cut's ribbon graph, the antibonding metric's
energy window) still builds on or was designed in explicit contrast to.

### Installation

```bash
pip install -r requirements.txt   # pymatgen
```

Nécessite Python ≥ 3.9 (utilise `pymatgen.io.lobster.outputs.Icohplist`).

### Format d'entrée attendu

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

### Utilisation

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

### Sortie

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

### Algorithme (résumé)

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

### Tests

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

### Exemples

`examples/dataset/` contient trois composés jouets illustrant les trois
cas ci-dessus (anisotrope, chemin indirect moins cher que la liaison
directe, direction déconnectée) :

```bash
python percolation_path.py --root examples/dataset --metric icohp --output /tmp/example_results.csv -v
```

### Analysis of `percolation_path.py` (missions #1-#2)

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

`analysis/REPORT_conventional_pilot.md` (mission #2): tested the pilot
report's hypothesis that the primitive/conventional cell choice explains
the near-zero correlation above, on the 6 pilot compounds. Verdict: the
choice does change the percolation weight substantially for compounds
whose cell actually differs (3x-55x, well past DFT noise), but the
strongest bond still never participates in the winning path, and the
weight gets *smaller* rather than larger as hypothesized — a bigger cell
just gives the minimum-weight search more long-range weak bonds to
exploit. Extension to the full dataset was not pursued based on this
result.

## Limites connues

- Suppose que `ICOHPLIST.lobster`/`ICOBILIST.lobster` proviennent du même
  calcul LOBSTER que la structure fournie (même ordre d'atomes, même
  maille) ; aucune vérification croisée de cohérence structure/POSCAR
  LOBSTER au-delà de la validation des indices d'atome.
- La complexité du Dijkstra par direction est `O(atomes)` exécutions sur un
  espace d'état de taille `atomes × (2·coord_bound+1)³` ; adapté à des
  mailles de quelques dizaines à quelques centaines d'atomes par composé,
  pas à des mailles géantes.
