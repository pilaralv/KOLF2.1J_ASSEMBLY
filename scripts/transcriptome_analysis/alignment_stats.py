from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import OUTDIR, PLOTDIR, STATDIR
from gff_utils import is_excluded, ref_key, sample_name


def parse_stats(path: Path) -> dict[str, float]:
    data: dict[str, float] = {}
    for line in path.read_text().splitlines():
        k, v = line.split("\t", 1)
        if k == "file":
            continue
        data[k] = float(v)
    return data


def write_alignment_summary(by_sample: dict) -> None:
    rows = []
    for sample in sorted(by_sample):
        refs = by_sample[sample]
        if "hg38" not in refs or "kolf" not in refs:
            continue
        hg, kf = refs["hg38"], refs["kolf"]
        rows.append(
            {
                "sample": sample,
                "hg38_error_rate": hg["error_rate"],
                "kolf_error_rate": kf["error_rate"],
                "delta_error_rate": kf["error_rate"] - hg["error_rate"],
                "hg38_supplementary": int(hg["supplementary"]),
                "kolf_supplementary": int(kf["supplementary"]),
                "delta_supplementary": int(kf["supplementary"] - hg["supplementary"]),
                "hg38_supp_per_primary": hg["supplementary"] / hg["primary"] if hg["primary"] else 0,
                "kolf_supp_per_primary": kf["supplementary"] / kf["primary"] if kf["primary"] else 0,
            }
        )

    if not rows:
        print("No complete hg38 vs KOLF alignment stat pairs yet.")
        return

    out = OUTDIR / "alignment_summary.tsv"
    cols = list(rows[0].keys())
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    samples = [r["sample"] for r in rows]
    x = np.arange(len(samples))
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    d_err = [r["delta_error_rate"] * 100 for r in rows]
    ax.bar(x, d_err, color=["#1b7837" if d < 0 else "#d73027" for d in d_err])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("delta error rate (KOLF - hg38, %)")
    ax.set_title("Alignment mismatch rate")

    ax = axes[1]
    d_supp = [r["delta_supplementary"] for r in rows]
    ax.bar(x, d_supp, color=["#1b7837" if d < 0 else "#d73027" for d in d_supp])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("delta supplementary alignments")
    ax.set_title("Alignment parsimony")

    fig.tight_layout()
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTDIR / "alignment_improvement_summary.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out} and alignment_improvement_summary.png ({len(rows)} pairs)")


def summarize_alignment() -> None:
    by_sample: dict[str, dict[str, dict[str, float]]] = {}
    for p in sorted(STATDIR.glob("*.stats.tsv")):
        bam = p.stem
        try:
            sample = sample_name(bam + ".bam")
            key = ref_key(bam)
        except ValueError:
            continue  # unknown reference prefix (e.g. CHM13)
        if is_excluded(sample):
            continue
        by_sample.setdefault(sample, {})[key] = parse_stats(p)

    write_alignment_summary(by_sample)
