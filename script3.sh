#!/usr/bin/env bash
set -euo pipefail

# ── Paramètres ────────────────────────────────────────────────────────────────
VCF="gq__genbank_filtered_cohort.vcf.gz"
REF="trimmed_genomic.fna"
GTF="genomicchrom.gtf"
SCG="9.txt"

# ── Dossiers de sortie ────────────────────────────────────────────────────────
mkdir -p vcfs consensus_fa cds_out logs

# ── Index de la référence (si absent) ────────────────────────────────────────
samtools faidx "$REF" 2>/dev/null || true

# ── Boucle sur tous les samples ───────────────────────────────────────────────
for SAMPLE in $(bcftools query -l "$VCF"); do

  echo ""
  echo "=============================="
  echo " Traitement : $SAMPLE"
  echo "=============================="

  # 1. Extraire les SNPs de ce sample uniquement
  bcftools view -v snps -s "$SAMPLE" -Oz \
    -o "vcfs/${SAMPLE}.vcf.gz" "$VCF"
  tabix -p vcf "vcfs/${SAMPLE}.vcf.gz"

  # 2. Générer le FASTA consensus personnalisé
  bcftools consensus \
    -f "$REF" \
    -s "$SAMPLE" \
    "vcfs/${SAMPLE}.vcf.gz" \
    > "consensus_fa/genome_${SAMPLE}.fasta"

  samtools faidx "consensus_fa/genome_${SAMPLE}.fasta"

  # 3. Extraire les CDS SCG (v3)
  #    ⚠️  La redirection 2> doit être EN DERNIER, après tous les arguments
  python3 extract_scg_cds_v3.py \
    --gtf      "$GTF" \
    --fasta    "consensus_fa/genome_${SAMPLE}.fasta" \
    --scg      "$SCG" \
    --sample   "$SAMPLE" \
    --out      "cds_out/${SAMPLE}_SCG_CDS.fasta" \
    --report   "cds_out/${SAMPLE}_report.tsv" \
    --excluded "cds_out/${SAMPLE}_excluded.tsv" \
    2> "logs/${SAMPLE}.extract_scg.log"

  echo "OK : CDS extraites pour $SAMPLE"
  echo "     → cds_out/${SAMPLE}_SCG_CDS.fasta"

done

echo ""
echo "=============================="
echo " Pipeline terminé"
echo " Résultats dans : cds_out/"
echo "=============================="