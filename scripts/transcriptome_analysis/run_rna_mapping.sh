#!/bin/bash
# ONT long-read cDNA analysis.
# Requirements: bash, samtools, stringtie, python3 (numpy, matplotlib)
#
# Inputs (in WORKDIR unless noted):
#   hg38.gff3                              GENCODE hg38 annotation
#   {hg38,KOLF2.1Jv1.1}_<sample>.bam       minimap2 -ax splice long-read BAMs
#   KOLF PAT/MAT GFF3s                     in GENESET_DIR (merged here)
#   excluded_samples.txt (optional)        substring patterns to skip
# Usage:
#   bash run_rna_mapping.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$(pwd)}"
CPUS="${CPUS:-32}"
FORCE="${FORCE:-0}"
SKIP_SUMMARY="${SKIP_SUMMARY:-0}"

GENESET_DIR="${GENESET_DIR:-${WORKDIR}}"
HG38_GFF="${WORKDIR}/hg38.gff3"
KOLF_PAT_GFF="${KOLF_PAT_GFF:-${GENESET_DIR}/KOLF_PAT_noY.gff3}"
KOLF_MAT_GFF="${KOLF_MAT_GFF:-${GENESET_DIR}/KOLF_MAT.gff3}"
KOLF_GFF="${WORKDIR}/KOLF_diploid_noY.gff3"
KOLF_ID_MAP="${WORKDIR}/KOLF_id_map.tsv"
QUANT_ROOT="${WORKDIR}/stringtie_quant"
STAT_ROOT="${WORKDIR}/alignment_stats_individual"

# ONT cDNA libraries to process (edit to match your cohort, or provide
# SAMPLE_LIST=<file> with one sample per line).
SAMPLES=(
  ASTRO_IND12_cDNA_PAU23737
  ASTRO_IND12_cDNA_PAU24722
  ASTRO_IND12_cDNA_PAU24790
  ASTRO_IND12_cDNA_PAU24904
  ASTRO_NYSCF_cDNA_PAU25463
  IPS_cDNA_PAM31312
  MGL_cDNA_PAM72486
  NEURON_cDNA_PAM31307
  NGN2_cDNA_PAM74535
  OLIGO_NYSCF_cDNA_PAU24958
  OLIGO_NYSCF_cDNA_PAU25487
)
if [[ -n "${SAMPLE_LIST:-}" && -f "${SAMPLE_LIST}" ]]; then
  mapfile -t SAMPLES < <(grep -vE '^\s*(#|$)' "${SAMPLE_LIST}")
fi

mkdir -p "$QUANT_ROOT" "$STAT_ROOT" "${WORKDIR}/reference_comparison"

merge_kolf_gff() {
  for f in "$KOLF_PAT_GFF" "$KOLF_MAT_GFF"; do
    if [[ ! -f "$f" ]]; then
      echo "[gff] ERROR: missing KOLF annotation: $f" >&2
      exit 1
    fi
  done
  if [[ -s "$KOLF_GFF" && "$FORCE" != "1" ]]; then
    echo "[gff] skip (exists): $KOLF_GFF"
    return 0
  fi
  echo "[gff] merge PAT/MAT annotations -> $(basename "$KOLF_GFF")"
  {
    echo "##gff-version 3"
    awk -F'\t' 'NF >= 9 && $1 !~ /^#/' "$KOLF_PAT_GFF"
    awk -F'\t' 'NF >= 9 && $1 !~ /^#/' "$KOLF_MAT_GFF"
  } > "$KOLF_GFF"
}

build_kolf_id_map() {
  if [[ -s "$KOLF_ID_MAP" && "$FORCE" != "1" ]]; then
    echo "[idmap] skip (exists): $KOLF_ID_MAP"
    return 0
  fi
  echo "[idmap] extract CAT gene_id -> source Ensembl ID"
  {
    echo -e "cat_gene_id\tensembl_base\tgene_name"
    awk -F'\t' '$3 == "gene" {
      n = split($9, a, ";");
      gid = ""; ensg = ""; gname = "";
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^gene_id=/)          gid = substr(a[i], 9);
        else if (a[i] ~ /^source_gene=/) ensg = substr(a[i], 13);
        else if (a[i] ~ /^gene_name=/)   gname = substr(a[i], 11);
      }
      sub(/\.[0-9]+(_[0-9]+)?$/, "", ensg);
      if (gid != "") print gid "\t" ensg "\t" gname;
    }' "$KOLF_PAT_GFF" "$KOLF_MAT_GFF"
  } > "$KOLF_ID_MAP"
  echo "[idmap] wrote $KOLF_ID_MAP"
}

stat_field() {
  awk -F'\t' -v key="$1" '$2 == key { print $3; exit }' "$2"
}

# Per-BAM alignment QC via samtools stats (primary counts + error/supp rates).
run_alignment_stats() {
  local bam_name="$1"
  local bam_path="${WORKDIR}/${bam_name}"
  local out="${STAT_ROOT}/${bam_name%.bam}.stats.tsv"

  if [[ -s "$out" && "$FORCE" != "1" ]]; then
    echo "  [stats] skip (exists): $out"
    return 0
  fi
  if [[ ! -f "$bam_path" ]]; then
    echo "  [stats] skip (missing BAM): $bam_path" >&2
    return 0
  fi

  local tmp
  tmp="$(mktemp "${STAT_ROOT}/.stats_XXXXXX")"
  trap 'rm -f "$tmp"' RETURN

  echo "  [stats] samtools stats $bam_name"
  samtools stats "$bam_path" > "$tmp"
  {
    echo -e "file\t${bam_name}"
    echo -e "primary\t$(samtools view -c -F 0x904 "$bam_path")"
    echo -e "reads_mq0\t$(stat_field "reads MQ0:" "$tmp")"
    echo -e "reads_mapped\t$(stat_field "reads mapped:" "$tmp")"
    echo -e "mismatches\t$(stat_field "mismatches:" "$tmp")"
    echo -e "error_rate\t$(stat_field "error rate:" "$tmp")"
    echo -e "supplementary\t$(stat_field "supplementary alignments:" "$tmp")"
    echo -e "bases_mapped_cigar\t$(stat_field "bases mapped (cigar):" "$tmp")"
  } > "$out"
  echo "  [stats] wrote $out"
}

# Reference-guided long-read quantification.
run_stringtie() {
  local bam_name="$1"
  local gff_path="$2"
  local bam_path="${WORKDIR}/${bam_name}"
  local outdir="${QUANT_ROOT}/${bam_name%.bam}"

  if [[ -s "${outdir}/gene_abundances.tab" && "$FORCE" != "1" ]]; then
    echo "  [stringtie] skip (exists): ${outdir}/gene_abundances.tab"
    return 0
  fi
  if [[ ! -f "$bam_path" ]]; then
    echo "  [stringtie] skip (missing BAM): $bam_path" >&2
    return 0
  fi

  mkdir -p "$outdir"
  echo "  [stringtie] -e -L on $bam_name with $(basename "$gff_path")"
  stringtie -e -G "$gff_path" -o "${outdir}/transcripts.gtf" \
    -A "${outdir}/gene_abundances.tab" -L -p "$CPUS" "$bam_path"
  echo "  [stringtie] wrote ${outdir}/gene_abundances.tab"
}

process_sample() {
  local sample="$1"
  echo ""
  echo "================================================================"
  echo "[$(date)] SAMPLE: $sample"
  echo "================================================================"
  for ref in hg38 KOLF2.1Jv1.1; do
    local bam_name="${ref}_${sample}.bam"
    local gff_path="$HG38_GFF"
    [[ "$ref" == "KOLF2.1Jv1.1" ]] && gff_path="$KOLF_GFF"
    echo "[$(date)] Reference: $ref"
    run_alignment_stats "$bam_name"
    run_stringtie "$bam_name" "$gff_path"
  done
}

echo "[$(date)] RNA mapping analysis (hg38 + KOLF diploid)"
echo "  workdir: $WORKDIR"
echo "  samples: ${#SAMPLES[@]}"
echo "  cpus: $CPUS"

merge_kolf_gff
build_kolf_id_map

for sample in "${SAMPLES[@]}"; do
  [[ -z "$sample" ]] && continue
  process_sample "$sample"
done

if [[ "$SKIP_SUMMARY" != "1" ]]; then
  echo "[$(date)] Python summaries (hg38 vs KOLF)"
  RNA_MAPPING_WORKDIR="$WORKDIR" python3 "${SCRIPT_DIR}/analysis/run_analysis.py" "$WORKDIR"
fi

echo "[$(date)] Done."
