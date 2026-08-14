#!/usr/bin/env python3
"""
extract_scg_cds_v3.py
─────────────────────────────────────────────────────────────────────────────
Pipeline en 3 étapes explicites :

  ÉTAPE 1 — Filtre pré-concat (per-segment)
      Chaque segment CDS est vérifié INDIVIDUELLEMENT dans son propre frame
      (phase GTF du segment). Si un stop interne est trouvé → gène exclu
      immédiatement, AVANT toute concaténation.
      → Capture les ~241 gènes avec stops pré-existants.

  ÉTAPE 2 — Concaténation + GTF frame
      Pour les gènes propres, on concatène et on applique la phase du premier
      segment. On note si des stops apparaissent (= artefacts de jonction).

  ÉTAPE 3 — Recherche du meilleur frame
      Si le frame GTF produit des stops après concaténation, on teste les 3
      offsets (0, 1, 2). Le meilleur est retenu. Si toujours des stops → exclu.

  REPORTING
      Les 476 "post_stop" de la v1 sont maintenant décomposés en :
        • pre_concat_stop   : stops dans les segments individuels (~241)
        • junction_stop_fixed   : stops de jonction corrigés par frame search
        • junction_stop_residual: stops de jonction non corrigés
      Rapport TSV avec colonnes : gene | len | gtf_phase | frame_used |
                                  n_segs | pre_concat_bad_segs | status
─────────────────────────────────────────────────────────────────────────────
Usage :
  python extract_scg_cds_v3.py \
      --gtf    genome.gtf     \
      --fasta  genome.fa      \
      --scg    scg_list.txt   \
      --out    out.fasta      \
      --sample SAMPLE_NAME    \
      --report report.tsv     \
      --excluded excluded.tsv
"""

import argparse
import re
import subprocess
from collections import defaultdict

# ─── Codons stop ────────────────────────────────────────────────────────────
STOPS = {"TAA", "TAG", "TGA"}


# ─── Utilitaires séquences ──────────────────────────────────────────────────

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def clean_seq(seq: str) -> str:
    return re.sub(r"[^ACGTN]", "N", seq.upper())


def internal_stop_count(seq: str) -> int:
    """Stops internes = tous les triplets SAUF le dernier."""
    L = (len(seq) // 3) * 3
    seq = seq[:L]
    if L < 6:
        return 0
    return sum(1 for i in range(0, L - 3, 3) if seq[i:i+3] in STOPS)


def strip_terminal_stop(seq: str) -> str:
    if len(seq) >= 3 and seq[-3:] in STOPS:
        return seq[:-3]
    return seq


# ─── Accès FASTA ─────────────────────────────────────────────────────────────

def faidx_fetch(fasta: str, chrom: str, start: int, end: int):
    region = f"{chrom}:{start}-{end}"
    p = subprocess.run(
        ["samtools", "faidx", fasta, region],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if p.returncode != 0:
        return None
    return "".join(p.stdout.splitlines()[1:]).strip()


# ─── Parsing GTF ─────────────────────────────────────────────────────────────

tx_re     = re.compile(r'transcript_id "([^"]+)"')
suffix_re = re.compile(r'(g\d+\.t1)$')


def parse_gtf(gtf: str, scg: set) -> dict:
    cds = defaultdict(list)
    for line in open(gtf):
        if line.startswith("#"):
            continue
        p = line.strip().split("\t")
        if len(p) != 9 or p[2] != "CDS":
            continue
        m = tx_re.search(p[8])
        if not m:
            continue
        tx = m.group(1)
        m2 = suffix_re.search(tx)
        if not m2:
            continue
        gid = m2.group(1)
        if gid not in scg:
            continue
        phase = int(p[7]) if p[7] in ("0", "1", "2") else 0
        cds[gid].append((p[0], int(p[3]), int(p[4]), p[6], phase))
    return cds


# ─── Écriture FASTA ──────────────────────────────────────────────────────────

def write_fasta(out, header: str, seq: str):
    out.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        out.write(seq[i:i+60] + "\n")


# ─── Récupération des segments ───────────────────────────────────────────────

def fetch_segments(segs: list, fasta: str):
    """
    Retourne list[(seq_brute, phase)] ou None si faidx échoue.
    La séquence est déjà revcomp si brin -.
    """
    pieces = []
    for chrom, s, e, st, ph in segs:
        seq = faidx_fetch(fasta, chrom, s, e)
        if seq is None:
            return None
        if st == "-":
            seq = revcomp(seq)
        pieces.append((clean_seq(seq), ph))
    return pieces


# ─── ÉTAPE 1 : Filtre pré-concat (per-segment) ───────────────────────────────

def check_segments_individual(pieces_with_phase: list) -> tuple:
    """
    Vérifie chaque segment INDIVIDUELLEMENT dans son propre frame (phase GTF).

    Retourne (has_stop: bool, n_bad_segs: int, details: list[str])
    details = liste des raisons par segment bad (pour le rapport).
    """
    n_bad   = 0
    details = []

    for idx, (seq, ph) in enumerate(pieces_with_phase):
        # Applique la phase de CE segment (pas une phase cumulée)
        seg = seq[ph:] if ph > 0 and len(seq) > ph else seq
        seg = strip_terminal_stop(seg)
        stops = internal_stop_count(seg)
        if stops > 0:
            n_bad += 1
            details.append(f"seg{idx}_ph{ph}_{stops}stops")

    return n_bad > 0, n_bad, details


# ─── ÉTAPE 3 : Meilleur frame après concaténation ────────────────────────────

def reconstruct_best_frame(raw: str, gtf_phase: int) -> tuple:
    """
    Teste les offsets [gtf_phase, 0, 1, 2] et retourne le meilleur.
    Retourne (cds, offset_used, n_stops_gtf, n_stops_best)
    n_stops_gtf = stops avec le frame GTF (avant toute correction)
    n_stops_best = stops avec le meilleur frame trouvé
    """
    # Stops avec le frame GTF (diagnostic de référence)
    cds_gtf    = strip_terminal_stop(raw[gtf_phase:]) if gtf_phase < len(raw) else raw
    stops_gtf  = internal_stop_count(cds_gtf)

    offsets    = list(dict.fromkeys([gtf_phase, 0, 1, 2]))
    best_cds   = cds_gtf
    best_off   = gtf_phase
    best_stops = stops_gtf

    for offset in offsets:
        if offset >= len(raw):
            continue
        cds   = strip_terminal_stop(raw[offset:])
        stops = internal_stop_count(cds)
        if stops < best_stops:
            best_stops = stops
            best_cds   = cds
            best_off   = offset
        if stops == 0:
            break

    return best_cds, best_off, stops_gtf, best_stops


# ─── Programme principal ─────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf",      required=True)
    ap.add_argument("--fasta",    required=True)
    ap.add_argument("--scg",      required=True)
    ap.add_argument("--out",      required=True)
    ap.add_argument("--sample",   required=True)
    ap.add_argument("--report",   required=True)
    ap.add_argument("--excluded", required=True)
    args = ap.parse_args()

    scg         = set(l.strip() for l in open(args.scg) if l.strip())
    cds_by_gene = parse_gtf(args.gtf, scg)
    stats       = defaultdict(int)

    with open(args.out,      "w") as out,  \
         open(args.report,   "w") as rep,  \
         open(args.excluded, "w") as excl:

        rep.write(
            "gene\tlen\tn_segs\tgtf_phase\tframe_used\t"
            "stops_gtf\tstops_best\tstatus\n"
        )
        excl.write("gene\treason\tdetail\n")

        for gid, segs in sorted(cds_by_gene.items()):
            stats["total"] += 1

            # ── Vérification du brin ────────────────────────────────────────
            strands = set(s[3] for s in segs)
            if len(strands) > 1:
                stats["bad_strand"] += 1
                excl.write(f"{gid}\tstrand_mismatch\t-\n")
                continue

            strand = segs[0][3]

            # ── Tri des segments ────────────────────────────────────────────
            if strand == "+":
                segs = sorted(segs, key=lambda x: x[1])
            else:
                segs = sorted(segs, key=lambda x: -x[2])

            gtf_phase = segs[0][4]   # phase du PREMIER segment

            # ── Récupération des séquences ──────────────────────────────────
            pieces = fetch_segments(segs, args.fasta)
            if pieces is None:
                stats["faidx_fail"] += 1
                excl.write(f"{gid}\tfaidx_fail\t-\n")
                continue

            # ══ ÉTAPE 1 : Filtre pré-concat (per-segment) ══════════════════
            # Chaque segment est vérifié dans son propre frame (phase GTF du
            # segment). Si stop trouvé → exclu IMMÉDIATEMENT, avant concat.
            has_stop, n_bad_segs, seg_details = check_segments_individual(pieces)
            if has_stop:
                stats["pre_concat_stop"] += 1
                detail = ";".join(seg_details)
                excl.write(f"{gid}\tpre_concat_stop\t{detail}\n")
                continue

            # ══ ÉTAPE 2 : Concaténation ══════════════════════════════════════
            # À ce stade, aucun segment individuel n'a de stop.
            # On concatène les séquences BRUTES (sans appliquer de phase ici,
            # la phase sera gérée dans reconstruct_best_frame).
            raw_cds = "".join(seq for seq, ph in pieces)

            # ══ ÉTAPE 3 : Meilleur frame ══════════════════════════════════════
            cds, frame_used, stops_gtf, stops_best = reconstruct_best_frame(
                raw_cds, gtf_phase
            )

            n_segs = len(segs)

            if stops_best == 0:
                # ── CDS propre → on l'écrit ──────────────────────────────────
                stats["kept"] += 1

                if frame_used == gtf_phase and stops_gtf == 0:
                    status = "gtf_ok"
                elif frame_used != gtf_phase:
                    status = f"frame_corrected_to_{frame_used}"
                    stats["frame_corrected"] += 1
                    if stops_gtf > 0:
                        stats["junction_stop_fixed"] += 1
                else:
                    status = "gtf_ok"

                write_fasta(out, f"{args.sample}|{gid}", cds)
                rep.write(
                    f"{gid}\t{len(cds)}\t{n_segs}\t{gtf_phase}\t{frame_used}\t"
                    f"{stops_gtf}\t0\t{status}\n"
                )

            else:
                # ── Stops résiduels malgré tous les frames → exclu ───────────
                stats["residual_stop"] += 1
                stats["junction_stop_residual"] += 1
                excl.write(
                    f"{gid}\tresidual_stop\t"
                    f"best_frame={frame_used}_stops={stops_best}_gtf_stops={stops_gtf}\n"
                )
                rep.write(
                    f"{gid}\t{len(cds)}\t{n_segs}\t{gtf_phase}\t{frame_used}\t"
                    f"{stops_gtf}\t{stops_best}\texcluded_residual\n"
                )

    # ── Résumé console ──────────────────────────────────────────────────────
    total_excl = (stats["bad_strand"] + stats["faidx_fail"]
                  + stats["pre_concat_stop"] + stats["residual_stop"])

    print("\n" + "═" * 55)
    print(f"  STATS — {args.sample}")
    print("═" * 55)
    print(f"  SCG dans le GTF              : {stats['total']}")
    print()
    print(f"  ── EXCLUSIONS ──────────────────────────────────────")
    print(f"  brin mixte                   : {stats['bad_strand']}")
    print(f"  faidx fail                   : {stats['faidx_fail']}")
    print(f"  [ÉTAPE 1] stop pré-concat    : {stats['pre_concat_stop']}"
          f"   ← segments individuels avec stops")
    print(f"  [ÉTAPE 3] stop résiduel      : {stats['residual_stop']}"
          f"   ← stops de jonction non corrigés")
    print(f"  total exclus                 : {total_excl}")
    print()
    print(f"  ── CONSERVÉS ───────────────────────────────────────")
    print(f"  TOTAL CONSERVÉS              : {stats['kept']}")
    print(f"    dont frame GTF OK          : {stats['kept'] - stats['frame_corrected']}")
    print(f"    dont frame corrigé         : {stats['frame_corrected']}"
          f"   (stops de jonction réparés)")
    print("═" * 55)
    print()
    print("  ── DÉCOMPOSITION des post_stop (v1) ──────────")
    print(f"  pre_concat_stop              : {stats['pre_concat_stop']}"
          f"   (attendu ~241)")
    print(f"  junction_stop_fixed          : {stats['junction_stop_fixed']}"
          f"   (artefacts de jonction corrigés)")
    print(f"  junction_stop_residual       : {stats['junction_stop_residual']}"
          f"   (artefacts non corrigés)")
    print(f"  total stops détectés         : "
          f"{stats['pre_concat_stop'] + stats['junction_stop_fixed'] + stats['junction_stop_residual']}")
  
    print("═" * 55 + "\n")


if __name__ == "__main__":
    main()