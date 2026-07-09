#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from pathlib import Path

if len(sys.argv) > 1:
    os.environ["RNA_MAPPING_WORKDIR"] = sys.argv[1]

sys.path.insert(0, str(Path(__file__).resolve().parent))

from alignment_stats import summarize_alignment
from tpm_comparison import run_tpm_comparison, summarize_tpm


def main() -> None:
    run_tpm_comparison()
    summarize_tpm()
    summarize_alignment()


if __name__ == "__main__":
    main()
