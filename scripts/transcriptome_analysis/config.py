from __future__ import annotations

import os
from pathlib import Path

WORKDIR = Path(os.environ.get("RNA_MAPPING_WORKDIR", ".")).resolve()

HG38_GFF = WORKDIR / "hg38.gff3"
KOLF_ID_MAP = WORKDIR / "KOLF_id_map.tsv"
QUANT = WORKDIR / "stringtie_quant"
STATDIR = WORKDIR / "alignment_stats_individual"
EXCLUDE_FILE = WORKDIR / "excluded_samples.txt"

OUTDIR = WORKDIR / "reference_comparison"
PLOTDIR = OUTDIR / "plots"
PAIRDIR = OUTDIR / "per_sample"
SUMMARY_IN = OUTDIR / "summary.tsv"
SUMMARY_OUT = OUTDIR / "tpm_summary.tsv"

REF_PREFIXES = {"hg38": "hg38_", "kolf": "KOLF2.1Jv1.1_"}

TPM_THRESHOLDS = (0.0, 1.0, 5.0, 10.0)
OUTLIER_MIN_TPM = 10.0
OUTLIER_MAD_Z = 4.0
OUTLIER_MIN_ABS_LOG2FC = 0.5
OUTLIER_LABEL_TOP = 12
