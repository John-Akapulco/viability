# Vérification indépendante des polymorphes CaO et TiO2 du manuscrit

Note technique — campagne du 2026-09-21/22, en appui au dossier de réponse aux
reviewers de "A Heuristic Approach to Rationalize Metastability of Materials
Above the Convex Hull by a Chemical-Bonding Criterion" (Reitz & Dronskowski,
ic-2026-04181q).

Objectif : recalculer indépendamment ΔE et ΔICOHP pour les trois
transformations du manuscrit concernant les composés ajoutés en révision
(eq. 5 : CaO sphalérite→rocksalt, eq. 6 : CaO CsCl→rocksalt, eq. 9 : TiO2
Pnma→rutile) avec une méthode volontairement alignée sur celle du
manuscrit, pour tester la robustesse du signe exobondic exigée par reviewer 2
("determine whether the sign remains robust to the choice of reference
state"). Les polymorphes CaO sont traités en PBEsol+D3(BJ), la fonctionnelle
du manuscrit pour les composés principaux du groupe ; TiO2 est traité
séparément en r2SCAN (dernière section ci-dessous), la fonctionnelle que le
manuscrit utilise spécifiquement pour ce composé.

## Méthode

- Fonctionnelle : PBEsol (`GGA = PS`) + D3 à amortissement Becke-Johnson
  (`IVDW = 12`) — choix motivé par la référence [56] du manuscrit (Grimme,
  Ehrlich, Goerigk, *J. Comput. Chem.* 2011, article sur la fonction
  d'amortissement, cohérent avec un choix BJ).
- Relaxation : `IBRION = 2`, `ISIF = 7` (volume seul ; les trois structures
  ont leurs positions atomiques fixées par symétrie — aucun degré de liberté
  interne à relaxer), `EDIFFG = -0.005`, `ISYM = 2` pendant la relaxation.
- Calcul de production (statique, compatible LOBSTER) : `IBRION = -1`,
  `NSW = 0`, `ISYM = -1`, `ENCUT = 700 eV` (valeur donnée dans la section
  Computational Details du manuscrit).
- Maillages k repris du Tableau S1 du SI : 13×13×13 (rocksalt, Fm-3m),
  11×11×11 (sphalérite, F-43m), 15×15×15 (CsCl, Pm-3m).
- POTCAR : `Ca_sv` (06Sep2000) + `O` (08Apr2002), identiques à ceux du
  calcul original du manuscrit (vérifié).
- Base LOBSTER : `Ca 3p 3s 4s` / `O 2p 2s` (pbeVASPfit2015), identique dans
  tous les calculs.
- Structures de départ : CsCl et sphalérite construites comme prototypes
  idéaux (positions Wyckoff fixes) à partir des distances Ca–O citées dans
  le texte du manuscrit ; rocksalt à partir de `mp-2605` (Materials
  Project). Les trois ont ensuite été relaxées en volume avec la méthode
  ci-dessus.
- Répertoires de calcul : `mp_dataset/structures/extension_CaOCsCl_lit`
  (committé), `verify_CaO_rocksalt`, `verify_CaO_sphalerite` (+
  `_relax_CaOCsCl` pour l'étape de relaxation) -- ces trois derniers sont
  des calculs de vérification ponctuels, non committés, gardés localement.
  Jobs SLURM : 61020-61028.

## Découverte préalable (avant cette campagne) : bug de propagation dans le manuscrit

Indépendamment de cette campagne, la comparaison réponse-aux-reviewers /
manuscrit surligné / SI a révélé que l'équation (5) et le Tableau S2 du SI
affichent **ΔICOHP = -75 kJ/mol** pour sphalérite→rocksalt, alors que la
section "Calculus of Bonding Energetics" du manuscrit, elle-même corrigée en
réponse à reviewer 2 (ICOHP rocksalt = -0.713 eV/Ca, et non -0.706 eV/Ca),
redonne en interne **-79.2 kJ/mol**. La correction du -0.706→-0.713 eV/Ca n'a
donc pas été répercutée dans l'éq. (5) ni dans la Table S2 — à corriger avant
resoumission.

## Découverte n°1 : le pipeline `extension_` du projet n'est pas en PBEsol+D3

Le pipeline utilisé pour l'ensemble du dataset (~600 composés,
`mp_dataset/prepare_vasp_lobster.py`, `INCAR_SETTINGS`) tourne en **PBE nu,
sans correction D3** (confirmé par OUTCAR : `GGA = --`, aucune ligne
`IVDW`), et non en PBEsol+D3 comme le manuscrit et comme certaines notes de
métadonnées locales l'affirmaient à tort (ex. `extension_TiO2Pnma_mp-754769/
mp_metadata.json`, note corrigée depuis). Ceci concerne uniquement le
pipeline `extension_`/dataset du projet — **le calcul original du
manuscrit, lui, est bien en PBEsol+D3** (vérifié directement dans l'OUTCAR de
`manuscript_CaO_rocksalt` : `GGA = Ps`, `IVDW = 12`).

## Découverte n°2 : paramètre de maille du rocksalt-CaO, écart avec le Tableau S1 du SI

| Polymorphe | Notre relaxation (PBEsol+D3-BJ) | SI Table S1 | Écart |
|---|---|---|---|
| CsCl (Pm-3m) | a = 2.86010 Å | a = 2.86011 Å | ~0% |
| Sphalérite (F-43m) | a = 5.14706 Å | a = 5.14699 Å | ~0% |
| **Rocksalt (Fm-3m)** | **a = 4.73037 Å** | **a = 4.77501 Å** | **-0.94%** |

Pour CsCl et sphalérite, la relaxation indépendante reproduit le Tableau S1
à la 4e-5e décimale. Pour le rocksalt, notre valeur (4.73037 Å) ne colle
**pas** au Tableau S1 (4.77501 Å), mais correspond quasi exactement à la
géométrie extraite directement du `vasprun.xml` du calcul original du
manuscrit conservé localement (`manuscript_CaO_rocksalt`, a ≈ 4.7305 Å,
distance Ca–O = 2.3652 Å contre les "2.39 Å" cités dans le texte du
manuscrit). Deux calculs indépendants (celui du manuscrit et le nôtre)
convergent donc vers la même valeur, différente de celle publiée dans le
Tableau S1 → **erreur probable dans le Tableau S1 du SI** (paramètre de
maille du rocksalt-CaO), à vérifier par les auteurs avant resoumission.

## Résultats finaux (méthode uniforme, structures fraîchement relaxées)

ICOHP par liaison Ca–O et totaux (plus proches voisins uniquement) :

| Polymorphe | ICOHP/liaison | Nb liaisons | ICOHP total |
|---|---|---|---|
| Rocksalt | -0.7370 eV | 6 | -426.7 kJ/mol |
| Sphalérite | -0.8648 eV | 4 | -333.8 kJ/mol |
| CsCl | -0.4339 eV | 8 | -334.9 kJ/mol |

**eq. (5) sphalérite → rocksalt**

| | Notre calcul | Manuscrit |
|---|---|---|
| ΔE | **-55.25 kJ/mol** | -55 kJ/mol (quasi exact) |
| ΔICOHP | -92.89 kJ/mol | -75 kJ/mol (périmé, Table S2/éq.5) / -79.2 kJ/mol (recalcul interne cohérent, section Calculus) |

**eq. (6) CsCl → rocksalt**

| | Notre calcul | Manuscrit |
|---|---|---|
| ΔE | **-82.16 kJ/mol** | -82 kJ/mol (quasi exact) |
| ΔICOHP | -91.78 kJ/mol | -116 kJ/mol |

## Interprétation

1. **Les énergies totales (ΔE) sont reproduites presque exactement** dans
   les deux réactions (écart <0.5%), ce qui valide fortement la géométrie et
   la fonctionnelle utilisées ici.
2. **Le signe exobondic est confirmé indépendamment** pour les trois
   polymorphes, à la fois par ΔE et par ΔICOHP (tous négatifs dans tous les
   cas) — conclusion qualitative du manuscrit solide, y compris pour les
   deux nouveaux polymorphes ajoutés en révision (CsCl-CaO, TiO2 Pnma).
3. **La magnitude du ΔICOHP reste décalée**, et le sens de l'écart
   s'inverse entre les deux réactions (notre valeur est plus grande en
   magnitude pour eq.5, plus petite pour eq.6) — le rocksalt, commun aux
   deux réactions comme produit, en est la source la plus probable : son
   ICOHP/liaison varie de -0.713 eV (texte du manuscrit) à -0.719 eV
   (relance LOBSTER sur les fichiers originaux du manuscrit) à -0.737 eV
   (notre relaxation indépendante) selon la version du calcul considérée,
   soit un spread de ~3.4% sans que géométrie, fonctionnelle, POTCAR, ENCUT
   ou maillage k n'en soient la cause identifiée (tous vérifiés identiques
   ou explicitement testés et écartés, y compris l'ENCUT qui a été testé à
   600 et 700 eV sans amélioration).
4. Piste non résolue : le résiduum semble se situer au niveau du calcul
   LOBSTER/COHP lui-même (convergence SCF résiduelle, fenêtre d'intégration
   COHP, ou une autre finesse numérique non identifiée) plutôt qu'au niveau
   DFT (VASP) — pas creusé davantage faute d'instruction spécifique.

## Recommandations pour la révision

- Corriger le paramètre de maille du rocksalt-CaO dans le Tableau S1 du SI
  (4.77501 Å semble erroné ; ~4.730 Å est ce que produit le calcul réel,
  confirmé indépendamment).
- Corriger l'éq. (5)/Table S2 : ΔICOHP(sphalérite→rocksalt) = -79.2 kJ/mol
  (pas -75 kJ/mol), pour rester cohérent avec la correction -0.706→-0.713
  eV/Ca déjà appliquée dans le texte suite à reviewer 2.
- Ne pas sur-interpréter la précision numérique du ΔICOHP à 2-3 chiffres
  significatifs : nos vérifications indépendantes montrent une
  reproductibilité de l'ordre de ±15-25 kJ/mol sur ces petites mailles
  ioniques très symétriques, même à géométrie/fonctionnelle/ENCUT/maillage
  identiques. Le signe, lui, est robuste dans tous les tests effectués.

## TiO2 Pnma → rutile (eq. 9) : vérification en r2SCAN

Contrairement à CaO, le manuscrit calcule spécifiquement TiO2 avec la
fonctionnelle méta-GGA r2SCAN (et non PBEsol+D3), parce que la structure
Pnma provient de Materials Project où elle a été relaxée en r2SCAN
(Computational Details du manuscrit). Une première passe de cette
vérification (2026-09-21) avait utilisé par erreur le PBE nu du pipeline
`extension_` du projet (voir Découverte n°1 ci-dessus), donnant un accord
correct sur le signe mais pas sur la fonctionnelle. Reprise en r2SCAN
propre (`METAGGA = R2scan`, `ENCUT = 700`, relaxation `IBRION = 2`/`ISIF =
3` avec `SYMPREC = 1e-4` pour lever une erreur de symétrie VASP sur le
rutile) :

| | SI Table S1 | Notre relaxation r2SCAN | Écart |
|---|---|---|---|
| Pnma (a, b, c) | 18.76727 / 2.96575 / 4.68009 Å | 18.76552 / 4.68258 / 2.96534 Å (axes b/c permutés) | ~0.01-0.05% |
| Rutile (a, c) | 4.60013 / 2.96018 Å | 4.60106 / 2.96016 Å | ~0.02% / ~0% |

| Réaction Pnma→rutile | Notre calcul (r2SCAN) | Manuscrit |
|---|---|---|
| ΔE | **-7.70 kJ/mol** | -8 kJ/mol (r2SCAN, quasi exact) / -10 kJ/mol (Materials Project) |
| ΔICOHP | **-9.48 kJ/mol** | -9 kJ/mol (quasi exact) |

**Accord essentiellement parfait**, à la fois sur ΔE et sur ΔICOHP, une fois
la bonne fonctionnelle (r2SCAN) utilisée — contrairement à CaO, ici aucun
résidu de magnitude significatif ne subsiste. Ceci confirme que le résidu
ICOHP observé sur les polymorphes CaO (section précédente) n'est pas un
artefact générique de notre pipeline LOBSTER, mais bien quelque chose de
spécifique à ce système (probablement lié au rocksalt-CaO, cf.
Interprétation ci-dessus) : quand fonctionnelle, géométrie, ENCUT et
maillage sont correctement alignés sur le manuscrit, l'accord peut être
quasi parfait, comme le montre TiO2.
