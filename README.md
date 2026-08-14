# Sélection positive sur les gènes à copie unique entre *Gonioctena quinquepunctata* et *G. intermedia*

Scripts développés pour le mémoire de Master en Bioinformatique et Modélisation (ULB) :
détection de la sélection positive entre deux espèces sœurs de chrysomèles, par le rapport
**Ka/Ks** et le test de **McDonald-Kreitman**, appliqués aux gènes à copie unique (orthologues 1:1).

**Auteur :** Bryan Derbrée — **Promoteur :** Patrick Mardulyn
(Laboratoire d'Évolution biologique et Écologie, ULB) — 2025-2026.

---

## Contenu du dépôt

| Fichier | Rôle |
|---|---|
| `script3.sh` | Script chef d'orchestre : pour chaque individu, extrait les SNP, génère le **consensus** (`bcftools`) et lance l'extraction des CDS. |
| `extract_scg_cds_v3.py` | Extrait les **CDS des gènes à copie unique** de chaque consensus, gère la phase de lecture et écarte les codons stop internes. |
| `prepare_mk.py` | **Regroupe** les CDS par gène (une séquence par individu, étiquetée par espèce) et applique un contrôle qualité. *Aucun réalignement* : les séquences sont pré-alignées par construction (pas d'indels). |
| `kaks_v2.py` | Calcule le **Ka/Ks** (KaKs_Calculator 2.0, méthode **YN**) par gène, avec test exact de Fisher et correction de Benjamini-Hochberg. |
| `resultats_figures.py` | Régénère l'ensemble des **figures et tableaux** du chapitre Résultats à partir des deux fichiers `.tsv`. |

> Le test de McDonald-Kreitman lui-même est réalisé avec l'outil **`mkado`** (Dn/Ds/Pn/Ps, α, DoS).

---

## Ordre d'exécution du pipeline

```
Données amont (génome de référence + annotation, appel des variants GATK,
orthologues OrthoFinder) — réutilisées de Lukicheva & Mardulyn (2021)
        │
        ▼
1. script3.sh              → consensus par individu + extract_scg_cds_v3.py
        │                     (sortie : {individu}_SCG_CDS.fasta)
        ▼
2. prepare_mk.py           → regroupement par gène
        │                     (sortie : mk_alignments/{gene}.fasta)
        ├───────────────┐
        ▼               ▼
3a. kaks_v2.py      3b. mkado
    → kaks_results.tsv   → mk_results_mkado.tsv
        └───────┬───────┘
                ▼
4. resultats_figures.py    → figures (.png) + tableaux (.csv)
```

---

## Dépendances

- **Python ≥ 3.8** : `pandas`, `numpy`, `matplotlib`
- **Outils bio-informatiques** : `bwa`, `GATK4`, `samtools`, `bcftools`, `tabix`,
  `KaKs_Calculator 2.0`, `mkado`

---

## Utilisation (exemples)

```bash
# 1. Consensus + extraction des CDS (tous les individus)
bash script3.sh

# 2. Regroupement par gène
python3 prepare_mk.py \
    --cds-dir cds_out/ \
    --samples samples.tsv \
    --out-dir mk_alignments/

# 3a. Ka/Ks (méthode YN)
python3 kaks_v2.py \
    --fasta-dir mk_alignments/ \
    --out kaks_results.tsv \
    --method YN

# 3b. Test de McDonald-Kreitman
#     (mkado sur les alignements par gène → mk_results_mkado.tsv)

# 4. Figures et tableaux
python3 resultats_figures.py kaks_results.tsv mk_results_mkado.tsv figures/
```

Le fichier `samples.tsv` associe chaque individu à son espèce (deux colonnes : `sample` et `population`).

---

## Données

Les données brutes (reséquençage de 20 individus, génome de référence de *G. quinquepunctata*)
proviennent de Lukicheva & Mardulyn (2021) et de Lukicheva, Flot & Mardulyn (2021).

## Licence

Code distribué sous licence **MIT** (voir le fichier `LICENSE`).
