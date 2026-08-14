#!/usr/bin/env python3
"""
prepare_mk.py
═════════════
Regroupe par gène les CDS consensus produites par extract_scg_cds_v3.py
(un fichier FASTA par individu) en un FASTA par gène, prêt pour le test de
McDonald-Kreitman (mkado) et pour kaks_v2.py.

Comme seuls les SNP ont été conservés à l'appel des variants (pas d'indels)
et que tous les consensus sont générés sur le même génome de référence, les
séquences d'un même gène ont la MÊME LONGUEUR et sont homologues position par
position : elles sont pré-alignées par construction. AUCUN réalignement n'est
donc effectué — les séquences sont simplement regroupées, étiquetées par
population, puis contrôlées.

Pour chaque gène SCG :
  1. Collecte les séquences CDS de tous les individus (pop1 + pop2)
  2. Étiquette chaque séquence : >{sample}__POP1  ou  >{sample}__POP2
  3. Contrôle qualité :
       - bases valides uniquement (A, C, G, T, N)
       - longueur identique pour tous les individus du gène
         (garantie d'homologie qui remplace l'alignement)
       - longueur multiple de trois (cadre intact)
       - longueur minimale (--min-codons)
       - nombre minimal d'individus par population (--min-seqs)
  4. Écrit le FASTA regroupé dans --out-dir/{gene_id}.fasta

Format des headers de sortie :
  >{sample}__POP1   ou   >{sample}__POP2
  → lisible directement par mkado, kaks_v2.py, egglib ou PopGenome

Usage :
  python3 prepare_mk.py \\
      --cds-dir  cds_out/          \\
      --samples  samples.tsv       \\
      --out-dir  mk_alignments/    \\
      [--min-seqs 2]               \\
      [--min-codons 1]

  # Format samples.tsv (tabulation) :
  #   sample          population
  #   INDIV_XA32      espece_A
  #   SAMPLE_BR91     espece_B
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict


# ─── Lecture FASTA ────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = {}
    current = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current is not None:
                    seqs[current] = ''.join(buf).upper()
                current = line[1:]
                buf = []
            else:
                buf.append(line)
    if current is not None:
        seqs[current] = ''.join(buf).upper()
    return seqs


def extract_gene_id(header):
    m = re.search(r'(g\d+\.t1)', header)
    return m.group(1) if m else header.split()[0]


# ─── Traitement d'un gène (regroupement + contrôle qualité, sans réalignement) ─

def process_gene(gid, seqs_pop1, seqs_pop2, out_dir, min_seqs, min_codons):
    """
    seqs_pop1 / seqs_pop2 : {sample_name: cds_sequence}
    Retourne : (gid, n_pop1, n_pop2, status)
    """
    n_p1 = len(seqs_pop1)
    n_p2 = len(seqs_pop2)

    if n_p1 < min_seqs:
        return gid, n_p1, n_p2, f'skip_pop1_too_few({n_p1})'
    if n_p2 < min_seqs:
        return gid, n_p1, n_p2, f'skip_pop2_too_few({n_p2})'

    # Étiquetage par population
    records = {}
    for sample, seq in seqs_pop1.items():
        records[f"{sample}__POP1"] = seq
    for sample, seq in seqs_pop2.items():
        records[f"{sample}__POP2"] = seq

    # QC — bases valides uniquement
    for label, seq in records.items():
        if set(seq) - set("ACGTN"):
            return gid, n_p1, n_p2, 'invalid_bases'

    # QC — longueur identique pour tous les individus (pré-alignés par construction)
    lengths = set(len(s) for s in records.values())
    if len(lengths) != 1:
        return gid, n_p1, n_p2, 'invalid_codon_aln_lengths'
    aln_len = lengths.pop()

    # QC — longueur multiple de trois (cadre intact)
    if aln_len % 3 != 0:
        return gid, n_p1, n_p2, f'invalid_codon_aln_not_multiple3({aln_len})'

    # QC — longueur minimale
    if aln_len < min_codons * 3:
        return gid, n_p1, n_p2, f'too_short({aln_len})'

    # Écriture du FASTA regroupé (séquences inchangées, sans alignement)
    out_path = os.path.join(out_dir, f"{gid}.fasta")
    with open(out_path, 'w') as fh:
        for label, seq in records.items():
            fh.write(f">{label}\n{seq}\n")

    return gid, n_p1, n_p2, 'OK'


# ─── Programme principal ──────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Regroupe les CDS par gène pour le test MK (sans réalignement)"
    )
    ap.add_argument('--cds-dir',    required=True,
                    help='Dossier contenant les {sample}_SCG_CDS.fasta')
    ap.add_argument('--samples',    required=True,
                    help='Fichier TSV : colonnes "sample" et "population"')
    ap.add_argument('--out-dir',    required=True,
                    help='Dossier de sortie pour les FASTA regroupés (un par gène)')
    ap.add_argument('--min-seqs',   type=int, default=2,
                    help='Séquences minimales par population par gène (défaut: 2)')
    ap.add_argument('--min-codons', type=int, default=1,
                    help='Longueur minimale en codons (défaut: 1)')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Découverte des fichiers CDS (un par individu) ─────────────────────────
    fasta_files = sorted(glob.glob(os.path.join(args.cds_dir, '*_SCG_CDS.fasta')))
    if not fasta_files:
        sys.exit(f"[ERREUR] Aucun fichier *_SCG_CDS.fasta dans {args.cds_dir}")

    available = {}
    for path in fasta_files:
        sample = os.path.basename(path).replace('_SCG_CDS.fasta', '')
        available[sample] = path

    # ── Lecture des métadonnées (sample → population) ─────────────────────────
    pop_labels = {}
    with open(args.samples) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            if len(parts) < 2:
                print(f"[WARN] ligne {i+1} ignorée (format invalide) : {line!r}")
                continue
            if i == 0 and parts[0].lower() in ('sample', 'nom', 'name', 'id'):
                continue
            pop_labels[parts[0].strip()] = parts[1].strip()

    unique_pops = sorted(set(pop_labels.values()))
    if len(unique_pops) != 2:
        sys.exit(f"[ERREUR] Il faut exactement 2 populations, trouvé : {unique_pops}")
    pop_a, pop_b = unique_pops
    print(f"Populations : '{pop_a}'  vs  '{pop_b}'")

    pop1_files, pop2_files = {}, {}
    for sample, pop in pop_labels.items():
        if sample not in available:
            print(f"[WARN] '{sample}' listé mais pas de fichier CDS — ignoré")
            continue
        (pop1_files if pop == pop_a else pop2_files)[sample] = available[sample]

    print(f"Individus pop1 ({pop_a}) : {len(pop1_files)}")
    print(f"Individus pop2 ({pop_b}) : {len(pop2_files)}")
    if not pop1_files or not pop2_files:
        sys.exit("[ERREUR] Une population est vide — vérifiez --samples")

    # ── Indexation par gène (LE REGROUPEMENT) ─────────────────────────────────
    print("\nRegroupement des CDS par gène ...")
    gene_pop1 = defaultdict(dict)   # {gene_id: {sample: seq}}
    gene_pop2 = defaultdict(dict)

    for sample, path in pop1_files.items():
        for header, seq in read_fasta(path).items():
            gene_pop1[extract_gene_id(header)][sample] = seq
    for sample, path in pop2_files.items():
        for header, seq in read_fasta(path).items():
            gene_pop2[extract_gene_id(header)][sample] = seq

    shared = sorted(set(gene_pop1) & set(gene_pop2))
    print(f"Gènes présents dans les 2 populations : {len(shared)}")
    print(f"Sortie : {args.out_dir}/\n")

    # ── Traitement gène par gène ──────────────────────────────────────────────
    results = []
    for k, gid in enumerate(shared, 1):
        results.append(
            process_gene(gid, gene_pop1[gid], gene_pop2[gid],
                         args.out_dir, args.min_seqs, args.min_codons)
        )
        if k % 500 == 0 or k == len(shared):
            print(f"  {k}/{len(shared)} gènes traités ...", flush=True)

    ok = sum(1 for r in results if r[3] == 'OK')
    skipped = [r for r in results if r[3] != 'OK']

    print()
    print("═" * 55)
    print("  RÉSUMÉ REGROUPEMENT (sans réalignement)")
    print("═" * 55)
    print(f"  Gènes traités       : {len(results)}")
    print(f"  FASTA écrits        : {ok}  →  {args.out_dir}/")
    if skipped:
        print(f"  Exclusions          : {len(skipped)}")
        for reason, count in Counter(r[3] for r in skipped).most_common():
            print(f"    {reason:<40} : {count}")
    print("═" * 55)


if __name__ == '__main__':
    main()
