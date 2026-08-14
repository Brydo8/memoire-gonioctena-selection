#!/usr/bin/env python3
"""
kaks_v2.py -- Pipeline Ka/Ks avec KaKs_Calculator 2.0
=======================================================
Format d'entree attendu :
  Un fichier FASTA par gene dans un meme dossier.
  Chaque fichier contient toutes les sequences (POP1 + POP2),
  avec le groupe encode dans le header apres "__" :
    >A00387_98_GW190421_AAGACCGT__POP1
    ATGCCC...
    >A00387_98_GW190421_AGCCTATC__POP2
    ATGCCC...
  Le nom du fichier = identifiant du gene (ex: g1234.t1.fasta -> g1234.t1).

Le fichier --samples est OPTIONNEL.
  Sans --samples : les groupes sont extraits automatiquement du suffixe "__GROUP".
  Avec --samples  : fichier TSV (sample<TAB>groupe) qui remplace la detection auto.
                    Le "sample" doit correspondre au premier mot du header FASTA.

Usage :
  python3 kaks_v2.py \
      --fasta-dir  /chemin/vers/fastas/ \
      --out        kaks_results.tsv \
      [--samples   samples.tsv] \
      [--method    YN] \
      [--kaks-bin  KaKs_Calculator] \
      [--outdir    ./kaks_workdir/] \
      [--min-len   50] \
      [--threads   4] \
      [--genetic-code 1]
=======================================================
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# -- Utilitaires --------------------------------------------------------------

def fmt(v, dec=4):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return 'NA'
    if isinstance(v, float):
        return f'{v:.{dec}f}'
    return str(v)

def bh_correction(pvalues):
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1] if not math.isnan(x[1]) else 2.0)
    adj = [float('nan')] * n
    prev = 1.0
    for rank, (orig, pval) in enumerate(reversed(indexed), 1):
        if math.isnan(pval):
            adj[orig] = float('nan')
            continue
        a = min(prev, pval * n / (n - rank + 1))
        adj[orig] = min(a, 1.0)
        prev = adj[orig]
    return adj

# -- Lecture FASTA ------------------------------------------------------------

def read_fasta(path):
    """Retourne {header: sequence} pour un fichier FASTA."""
    seqs = {}
    current_id = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(chunks).upper()
                current_id = line[1:]   # Conserver le header COMPLET (sans le >)
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(chunks).upper()
    return seqs

def extract_group_from_header(header):
    """
    Extrait le groupe encode dans le header apres '__'.
    Ex: "A00387_98_GW190421_AAGACCGT__POP1" -> "POP1"
    Retourne None si aucun '__' trouve.
    """
    if '__' in header:
        return header.rsplit('__', 1)[1].strip()
    return None

# -- Fichier de correspondance (optionnel) ------------------------------------

def load_sample_map(path):
    """
    Lit le fichier TSV individu->groupe.
    Retourne {sample: groupe}.
    """
    sample_to_group = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            if len(parts) < 2:
                continue
            sample, group = parts[0].strip(), parts[1].strip()
            if sample.lower() in ('sample', 'individu', 'individual', 'id', 'name'):
                continue
            sample_to_group[sample] = group
    return sample_to_group

# -- Chargement (un fichier par gene) ----------------------------------------

def load_genes_by_file(fasta_dir, sample_to_group=None):
    """
    Charge les genes depuis un dossier (un fichier FASTA = un gene).
    Chaque FASTA contient les sequences de tous les individus (POP1 + POP2).

    Assignation du groupe :
      - Si sample_to_group fourni : le premier mot du header est cherche dans le dict.
      - Sinon : le groupe est extrait du suffixe "__GROUP" du header.

    Retourne :
        {gene_id: {group: [seq1, seq2, ...]}}
    """
    EXTENSIONS = {'.fasta', '.fa', '.fna', '.fas'}
    gene_data = {}
    n_unknown = 0
    fasta_dir = Path(fasta_dir)

    fasta_files = sorted(
        p for p in fasta_dir.iterdir()
        if p.suffix.lower() in EXTENSIONS
    )

    if not fasta_files:
        print(f'[ERREUR] Aucun fichier FASTA trouve dans {fasta_dir}', file=sys.stderr)
        sys.exit(1)

    for fasta_file in fasta_files:
        gene_id = fasta_file.stem      # g1234.t1.fasta -> "g1234.t1"
        seqs = read_fasta(fasta_file)
        if not seqs:
            continue

        groups = defaultdict(list)
        for header, seq in seqs.items():
            group = None

            if sample_to_group is not None:
                # Mode samples.tsv : match sur le premier mot du header
                sample = header.split()[0]
                group = sample_to_group.get(sample)
                if group is None:
                    for k, g in sample_to_group.items():
                        if k.lower() == sample.lower():
                            group = g
                            break
            else:
                # Mode auto : groupe apres "__" dans le header
                group = extract_group_from_header(header)

            if group is not None:
                groups[group].append(seq)
            else:
                n_unknown += 1

        if groups:
            gene_data[gene_id] = dict(groups)

    if n_unknown > 0:
        print(f'  [ATTENTION] {n_unknown} sequences sans groupe detecte.', file=sys.stderr)
        if sample_to_group is None:
            print('  -> Verifiez que les headers contiennent "__GROUP" (ex: ...BARCODE__POP1)',
                  file=sys.stderr)
    return gene_data

# -- Consensus ----------------------------------------------------------------

def majority_consensus(sequences):
    """Vote majoritaire position par position. Ignore gaps et N."""
    if not sequences:
        return None
    seqs = [s.replace(' ', '') for s in sequences if len(s) > 3]
    if not seqs:
        return None
    lengths = Counter(len(s) for s in seqs)
    modal_len = lengths.most_common(1)[0][0]
    valid = [s for s in seqs if abs(len(s) - modal_len) <= 3]
    if not valid:
        valid = seqs
    ref_len = min(len(s) for s in valid)
    ref_len = ref_len - (ref_len % 3)
    if ref_len < 3:
        return None
    valid = [s[:ref_len] for s in valid]
    consensus = []
    for i in range(ref_len):
        bases = [s[i] for s in valid if i < len(s) and s[i] not in ('-', 'N', 'n')]
        consensus.append(Counter(bases).most_common(1)[0][0] if bases else 'N')
    return ''.join(consensus)

# -- Nettoyage CDS ------------------------------------------------------------

STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def clean_cds(seq, min_codons=50):
    seq = seq.upper().replace(' ', '').replace('-', '')
    seq = seq[:len(seq) - (len(seq) % 3)]
    if len(seq) < min_codons * 3:
        return None, f'trop_court ({len(seq)} nt)'
    if seq[-3:] in STOP_CODONS:
        seq = seq[:-3]
    for i in range(0, len(seq) - 3, 3):
        if seq[i:i+3] in STOP_CODONS:
            return None, f'stop_interne (pos {i})'
    if len(seq) < min_codons * 3:
        return None, f'trop_court_apres_nettoyage ({len(seq)} nt)'
    return seq, 'OK'

# -- Nettoyage codons ambigus (N/gaps) ----------------------------------------

def strip_ambiguous_codons(seq1, seq2):
    """
    Supprime les codons avec gaps ou N dans l'un ou l'autre consensus.
    Les sequences doivent avoir la meme longueur (pre-alignees).
    """
    assert len(seq1) == len(seq2), f'Longueurs differentes : {len(seq1)} vs {len(seq2)}'
    trim = len(seq1) - (len(seq1) % 3)
    seq1, seq2 = seq1[:trim], seq2[:trim]
    c1, c2 = [], []
    for i in range(0, len(seq1), 3):
        a, b = seq1[i:i+3], seq2[i:i+3]
        if '-' not in a and '-' not in b and 'N' not in a and 'N' not in b:
            c1.append(a)
            c2.append(b)
    return ''.join(c1), ''.join(c2)

# -- Format AXT ---------------------------------------------------------------

def write_axt(pairs, axt_path):
    with open(axt_path, 'w') as f:
        for gene_id, s1, s2 in pairs:
            f.write(f'{gene_id}\n{s1}\n{s2}\n\n')

# -- KaKs_Calculator ----------------------------------------------------------

def run_kaks_calculator(kaks_bin, axt_path, out_path, method='YN', genetic_code=1):
    cmd = [kaks_bin, '-i', str(axt_path), '-o', str(out_path),
           '-m', method, '-c', str(genetic_code)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f'  [ERREUR] stdout: {result.stdout[:200]}', file=sys.stderr)
            print(f'  [ERREUR] stderr: {result.stderr[:200]}', file=sys.stderr)
            return False
        return True
    except subprocess.TimeoutExpired:
        print('  [ERREUR] KaKs_Calculator timeout (> 1h)', file=sys.stderr)
        return False
    except FileNotFoundError:
        print(f'  [ERREUR] Binaire introuvable : {kaks_bin}', file=sys.stderr)
        return False

def parse_kaks_output(out_path):
    results = {}
    if not Path(out_path).exists():
        return results
    with open(out_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('Sequence', '').strip()
            if not gene:
                continue
            def sf(key):
                try:
                    v = float(row.get(key, 'NA') or 'NA')
                    return v if not math.isnan(v) else float('nan')
                except (ValueError, TypeError):
                    return float('nan')
            ka = sf('Ka')
            ks = sf('Ks')
            ratio = sf('Ka/Ks')
            if math.isnan(ratio) and not (math.isnan(ka) or math.isnan(ks)):
                ratio = ka / ks if ks > 0 else float('inf')
            results[gene] = {
                'Ka': ka, 'Ks': ks, 'Ka_Ks': ratio,
                'p_fisher': sf('P-Value(Fisher)'),
                'S_sites': sf('S-Sites'),
                'N_sites': sf('N-Sites'),
            }
    return results

# -- Traitement d'un gene (parallelisable) ------------------------------------

def process_gene(args):
    """
    Consensus par groupe (sequences pre-alignees) -> strip N/gaps -> clean CDS.
    Pas d'alignement MAFFT : les sequences sont deja alignees.
    """
    gene, seqs_grp1, seqs_grp2, min_codons = args

    c1 = majority_consensus(seqs_grp1)
    c2 = majority_consensus(seqs_grp2)
    if c1 is None or c2 is None:
        return gene, None, None, 'consensus_impossible'

    if len(c1) != len(c2):
        return gene, None, None, f'longueur_inegale_consensus ({len(c1)} vs {len(c2)})'

    c1, c2 = strip_ambiguous_codons(c1, c2)

    if len(c1) < min_codons * 3:
        return gene, None, None, f'trop_court ({len(c1)} nt)'

    c1, st1 = clean_cds(c1, min_codons)
    if c1 is None:
        return gene, None, None, st1
    c2, st2 = clean_cds(c2, min_codons)
    if c2 is None:
        return gene, None, None, st2

    min_len = min(len(c1), len(c2))
    min_len = min_len - (min_len % 3)
    c1, c2 = c1[:min_len], c2[:min_len]

    if len(c1) < min_codons * 3:
        return gene, None, None, f'trop_court_post_nettoyage ({len(c1)} nt)'

    return gene, c1, c2, 'OK'

# -- Programme principal ------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description='Pipeline Ka/Ks (un FASTA par gene, groupes dans les headers)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples de headers supportes :
  >A00387_98_GW190421_AAGACCGT__POP1   ->  groupe = POP1  (detection auto)
  >QUI_individu3                         ->  utiliser --samples pour assigner le groupe

Sans --samples, le groupe est extrait automatiquement du suffixe '__GROUP'.
        """
    )
    ap.add_argument('--fasta-dir',  required=True,
                    help='Dossier FASTA (un fichier = un gene)')
    ap.add_argument('--samples',    default=None,
                    help='[OPTIONNEL] Fichier TSV sample->groupe. '
                         'Si absent, le groupe est lu depuis "__GROUP" dans les headers.')
    ap.add_argument('--kaks-bin',   default='KaKs_Calculator',
                    help='Binaire KaKs_Calculator (defaut: dans le PATH)')
    ap.add_argument('--method',     default='YN',
                    help='Methode Ka/Ks : NG, YN, MYN, GY, MA, MS, GMYN... (defaut: YN)')
    ap.add_argument('--genetic-code', type=int, default=1,
                    help='Code genetique NCBI (defaut: 1=Standard, 2=Vert_Mito, 5=Invert_Mito)')
    ap.add_argument('--out',        default='kaks_results.tsv',
                    help='Fichier TSV de sortie')
    ap.add_argument('--outdir',     default='./kaks_workdir',
                    help='Dossier de travail temporaire')
    ap.add_argument('--min-len',    type=int, default=50,
                    help='Longueur minimale en codons (defaut: 50)')
    ap.add_argument('--threads',    type=int, default=4,
                    help='Threads pour les consensus (defaut: 4)')
    ap.add_argument('--keep-workdir', action='store_true',
                    help='Conserver le dossier de travail')
    args = ap.parse_args()

    # -- Fichier de correspondance (optionnel)
    sample_to_group = None
    if args.samples is not None:
        sample_to_group = load_sample_map(args.samples)
        print(f'  Correspondance : {args.samples} ({len(sample_to_group)} individus)')

    # -- Binaire KaKs_Calculator
    kaks_bin = shutil.which(args.kaks_bin) or (args.kaks_bin if os.path.isfile(args.kaks_bin) else None)
    if kaks_bin is None:
        print(f'[ERREUR] KaKs_Calculator introuvable : {args.kaks_bin}')
        print('  Verifiez que conda est active ou fournissez --kaks-bin /chemin/complet')
        sys.exit(1)

    workdir = Path(args.outdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print('================================================')
    print('  Pipeline Ka/Ks -- KaKs_Calculator 2.0')
    print('================================================')
    print(f'  Dossier FASTA  : {args.fasta_dir}')
    print(f'  Methode        : {args.method}  |  Code genetique : {args.genetic_code}')
    print(f'  KaKs_Calculator: {kaks_bin}')
    print(f'  Min codons     : {args.min_len}  |  Threads : {args.threads}')
    if sample_to_group is None:
        print('  Groupes        : detection auto depuis le suffixe "__GROUP" des headers')

    # -- [1/5] Chargement
    print('\n[1/5] Chargement des FASTA (un fichier = un gene)...')
    gene_data = load_genes_by_file(args.fasta_dir, sample_to_group)
    print(f'  Fichiers FASTA lus : {len(gene_data)}')

    # Detecter les groupes presents
    all_groups = set()
    for d in gene_data.values():
        all_groups.update(d.keys())

    if len(all_groups) < 2:
        print(f'[ERREUR] Moins de 2 groupes detectes : {all_groups}')
        print('  Verifiez que les headers contiennent "__GROUP" (ex: ...BARCODE__POP1)')
        sys.exit(1)

    group_names = sorted(all_groups)
    if len(group_names) > 2:
        print(f'[ATTENTION] Plus de 2 groupes : {group_names}. Seuls {group_names[:2]} utilises.')
        group_names = group_names[:2]

    grp1, grp2 = group_names

    # Filtrer genes avec les 2 groupes
    valid_genes = {
        g: d for g, d in gene_data.items()
        if grp1 in d and grp2 in d
    }
    n_missing = len(gene_data) - len(valid_genes)

    print(f'  Groupes detectes : {grp1} ({sum(1 for d in gene_data.values() if grp1 in d)} genes) '
          f'et {grp2} ({sum(1 for d in gene_data.values() if grp2 in d)} genes)')

    # Afficher le nombre d'individus dans le premier fichier
    first_gene = next(iter(valid_genes.values())) if valid_genes else {}
    n1 = len(first_gene.get(grp1, []))
    n2 = len(first_gene.get(grp2, []))
    print(f'  Individus / gene : {grp1}={n1}, {grp2}={n2} (premier gene)')
    print(f'  Genes avec les 2 groupes : {len(valid_genes)}')
    if n_missing:
        print(f'  Genes avec groupe(s) manquant(s) : {n_missing} (ignores)')

    if not valid_genes:
        print('\n[ERREUR] Aucun gene avec les deux groupes.')
        sys.exit(1)

    # -- [2/5] Consensus
    print(f'\n[2/5] Calcul des consensus ({args.threads} threads)...')
    job_args = [
        (gene, data[grp1], data[grp2], args.min_len)
        for gene, data in sorted(valid_genes.items())
    ]

    pairs_ok = []
    skipped  = {}
    n_done   = 0

    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_gene, a): a[0] for a in job_args}
        for fut in as_completed(futures):
            gene, s1, s2, status = fut.result()
            n_done += 1
            if n_done % 500 == 0:
                print(f'  {n_done}/{len(valid_genes)} genes traites...', flush=True)
            if status == 'OK':
                pairs_ok.append((gene, s1, s2))
            else:
                skipped[gene] = status

    print(f'  Paires pretes : {len(pairs_ok)}')
    print(f'  Genes exclus  : {len(skipped)}')
    for reason, n in Counter(skipped.values()).most_common():
        print(f'    {reason}: {n}')

    if not pairs_ok:
        print('\n[ERREUR] Aucune paire disponible pour Ka/Ks.')
        sys.exit(1)

    # -- [3/5] Fichier AXT
    print('\n[3/5] Ecriture du fichier AXT...')
    axt_path = workdir / 'pairs.axt'
    write_axt(pairs_ok, axt_path)
    print(f'  {axt_path} ({len(pairs_ok)} paires)')

    # -- [4/5] KaKs_Calculator
    print(f'\n[4/5] Lancement de KaKs_Calculator (methode {args.method})...')
    kaks_out = workdir / f'kaks_{args.method}.txt'
    if not run_kaks_calculator(kaks_bin, axt_path, kaks_out, args.method, args.genetic_code):
        print('[ERREUR] KaKs_Calculator a echoue.')
        sys.exit(1)
    print(f'  Resultats bruts : {kaks_out}')

    # -- [5/5] Parsing + BH
    print('\n[5/5] Parsing et correction BH...')
    kaks_results = parse_kaks_output(kaks_out)
    print(f'  Genes avec Ka/Ks : {len(kaks_results)}')

    final_rows = []
    for gene, s1, s2 in pairs_ok:
        r = kaks_results.get(gene, {})
        ka    = r.get('Ka',       float('nan'))
        ks    = r.get('Ks',       float('nan'))
        ratio = r.get('Ka_Ks',    float('nan'))
        pval  = r.get('p_fisher', float('nan'))
        if math.isnan(ratio) or math.isinf(ratio):
            status = 'indefini'
        elif ks == 0:
            status = 'Ks_zero'
        else:
            status = 'OK'
        final_rows.append({
            'gene': gene, 'n_codons': len(s1) // 3,
            'Ka': ka, 'Ks': ks, 'Ka_Ks': ratio, 'p_fisher': pval,
            'S_sites': r.get('S_sites', float('nan')),
            'N_sites': r.get('N_sites', float('nan')),
            'method': args.method, 'status': status,
        })

    for gene, reason in skipped.items():
        final_rows.append({
            'gene': gene, 'n_codons': 0,
            'Ka': float('nan'), 'Ks': float('nan'),
            'Ka_Ks': float('nan'), 'p_fisher': float('nan'),
            'S_sites': float('nan'), 'N_sites': float('nan'),
            'method': args.method, 'status': reason,
        })

    ok_rows = [r for r in final_rows if r['status'] == 'OK' and not math.isnan(r['p_fisher'])]
    adj = bh_correction([r['p_fisher'] for r in ok_rows])
    for r, padj in zip(ok_rows, adj):
        r['p_fisher_adj'] = padj
        ratio = r['Ka_Ks']
        if math.isnan(ratio):
            r['interpretation'] = 'indefini'
        elif ratio > 1:
            r['interpretation'] = 'positive_selection'
        elif ratio < 1:
            r['interpretation'] = 'purifying_selection'
        else:
            r['interpretation'] = 'neutral'
    for r in final_rows:
        if 'p_fisher_adj'   not in r: r['p_fisher_adj']   = float('nan')
        if 'interpretation' not in r: r['interpretation'] = 'NA'

    cols = ['gene', 'n_codons', 'S_sites', 'N_sites', 'Ka', 'Ks', 'Ka_Ks',
            'p_fisher', 'p_fisher_adj', 'interpretation', 'method', 'status']

    final_rows.sort(key=lambda r: (
        0 if r['status'] == 'OK' else 1,
        -(r['Ka_Ks'] if not math.isnan(r.get('Ka_Ks', float('nan'))) else 0)
    ))

    with open(args.out, 'w', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for r in final_rows:
            f.write('\t'.join(fmt(r.get(c)) for c in cols) + '\n')

    n_pos  = sum(1 for r in ok_rows if r['interpretation'] == 'positive_selection')
    n_pur  = sum(1 for r in ok_rows if r['interpretation'] == 'purifying_selection')
    n_neut = sum(1 for r in ok_rows if r['interpretation'] == 'neutral')
    n_sig  = sum(1 for r in ok_rows if not math.isnan(r['p_fisher_adj']) and r['p_fisher_adj'] < 0.05)

    print()
    print('================================================')
    print('  RESUME Ka/Ks')
    print('================================================')
    print(f'  {grp1} vs {grp2}')
    print(f'  Genes analyses            : {len(ok_rows)}')
    print(f'  Ka/Ks > 1 (sel. positive) : {n_pos}')
    print(f'  Ka/Ks < 1 (sel. purif.)   : {n_pur}')
    print(f'  Ka/Ks = 1 (neutre)        : {n_neut}')
    print(f'  Significatifs BH (q<0.05) : {n_sig}')
    print()

    if n_pos > 0:
        top = sorted(ok_rows, key=lambda r: r['Ka_Ks'] if not math.isnan(r['Ka_Ks']) else 0,
                     reverse=True)[:10]
        print('  Top 10 Ka/Ks :')
        print(f'  {"Gene":15s}  {"Ka/Ks":>7}  {"Ka":>7}  {"Ks":>7}  {"p_adj":>8}')
        for r in top:
            print(f'  {r["gene"]:15s}  {fmt(r["Ka_Ks"]):>7}  {fmt(r["Ka"]):>7}  '
                  f'{fmt(r["Ks"]):>7}  {fmt(r["p_fisher_adj"]):>8}')

    print(f'\n  Resultats complets : {args.out}')
    print(f'  Methode utilisee   : {args.method}')

    if not args.keep_workdir:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print(f'  Fichiers intermediaires : {workdir}')

    print()
    print('================================================')


if __name__ == '__main__':
    main()