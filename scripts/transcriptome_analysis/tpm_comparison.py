from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import (
    OUTDIR,
    OUTLIER_LABEL_TOP,
    OUTLIER_MAD_Z,
    OUTLIER_MIN_ABS_LOG2FC,
    OUTLIER_MIN_TPM,
    PAIRDIR,
    PLOTDIR,
    QUANT,
    SUMMARY_IN,
    SUMMARY_OUT,
    TPM_THRESHOLDS,
)
from gff_utils import (
    gene_base_id,
    is_excluded,
    is_named_protein_coding,
    load_excluded_patterns,
    load_hg38_gene_tpm,
    load_kolf_gene_tpm,
    load_kolf_id_map,
    load_named_protein_coding_genes,
    sample_name,
)


def pearson(xs: np.ndarray, ys: np.ndarray) -> float:
    if len(xs) < 2:
        return float("nan")
    return float(np.corrcoef(xs, ys)[0, 1])


def spearman(xs: np.ndarray, ys: np.ndarray) -> float:
    def ranks(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a)
        r = np.empty_like(a, dtype=float)
        i = 0
        while i < len(a):
            j = i
            while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
                j += 1
            avg = (i + j + 2) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    return pearson(ranks(xs), ranks(ys))


def count_genes_above(tpm: dict[str, dict], thresh: float, *, kolf_raw: bool = False) -> int:
    n = 0
    for g in tpm.values():
        val = float(g["tpm_raw"]) if kolf_raw else float(g.get("tpm", g.get("tpm_raw", 0.0)))
        if val > thresh:
            n += 1
    return n


def kolf_tpm_for_detection(kf: dict[str, dict], gene: str) -> float:
    return float(kf[gene]["tpm_raw"])


def gain_loss(genes: list[str], hg: dict[str, dict], kf: dict[str, dict]) -> tuple[int, int]:
    gained = lost = 0
    for g in genes:
        if kolf_tpm_for_detection(kf, g) > 1.0 and float(hg[g]["tpm"]) < 0.1:
            gained += 1
        if float(hg[g]["tpm"]) > 1.0 and kolf_tpm_for_detection(kf, g) < 0.1:
            lost += 1
    return gained, lost


def outlier_mask(diff_log: np.ndarray, expressed: np.ndarray) -> np.ndarray:
    mask = np.zeros(len(diff_log), dtype=bool)
    if not expressed.any():
        return mask
    d = diff_log[expressed]
    med = float(np.median(d))
    mad = float(np.median(np.abs(d - med)))
    if mad == 0:
        return mask
    scaled = 1.4826 * mad
    z = np.abs(diff_log - med) / scaled
    hit = (z > OUTLIER_MAD_Z) & (np.abs(diff_log) > OUTLIER_MIN_ABS_LOG2FC)
    mask[expressed] = hit[expressed]
    return mask


def _write_pair_tables(sample: str, pair_rows: list[dict]) -> None:
    PAIRDIR.mkdir(parents=True, exist_ok=True)

    with (PAIRDIR / f"{sample}_tpm_pairs.tsv").open("w") as fh:
        fh.write("gene_id\tgene_name\thg38_tpm\tkolf_tpm_raw\tdelta_tpm_raw\tlog2fc_raw\tis_outlier\n")
        for row in pair_rows:
            fh.write(
                f"{row['gene_id']}\t{row['gene_name']}\t"
                f"{row['hg38_tpm']:.6f}\t{row['kolf_tpm_raw']:.6f}\t"
                f"{row['delta_tpm_raw']:.6f}\t{row['log2fc_raw']:.6f}\t{row['is_outlier']}\n"
            )

    outlier_rows = sorted(
        (r for r in pair_rows if r["is_outlier"]),
        key=lambda r: abs(r["log2fc_raw"]),
        reverse=True,
    )
    with (PAIRDIR / f"{sample}_tpm_outliers.tsv").open("w") as fh:
        fh.write("gene_id\tgene_name\thg38_tpm\tkolf_tpm_raw\tdelta_tpm_raw\tlog2fc_raw\n")
        for row in outlier_rows:
            fh.write(
                f"{row['gene_id']}\t{row['gene_name']}\t"
                f"{row['hg38_tpm']:.6f}\t{row['kolf_tpm_raw']:.6f}\t"
                f"{row['delta_tpm_raw']:.6f}\t{row['log2fc_raw']:.6f}\n"
            )

    top = sorted(pair_rows, key=lambda r: r["abs_delta_tpm"], reverse=True)[:100]
    with (PAIRDIR / f"{sample}_top_delta_tpm.tsv").open("w") as fh:
        fh.write("gene_id\tgene_name\thg38_tpm\tkolf_tpm_raw\tdelta_tpm_raw\tis_outlier\n")
        for row in top:
            fh.write(
                f"{row['gene_id']}\t{row['gene_name']}\t"
                f"{row['hg38_tpm']:.6f}\t{row['kolf_tpm_raw']:.6f}\t"
                f"{row['delta_tpm_raw']:.6f}\t{row['is_outlier']}\n"
            )


def _plot_pair_scatter(
    sample: str,
    hg_log: np.ndarray,
    kf_log: np.ndarray,
    diff_log_all: np.ndarray,
    expressed: np.ndarray,
    outliers: np.ndarray,
    pair_rows: list[dict],
    md: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    mask = expressed
    inlier = mask & ~outliers
    ax = axes[0]
    if inlier.any():
        ax.scatter(hg_log[inlier], kf_log[inlier], s=1, alpha=0.12, c="#2166ac", rasterized=True)
    out = mask & outliers
    if out.any():
        ax.scatter(
            hg_log[out], kf_log[out], s=28, alpha=0.85, c="#ff7f00",
            edgecolors="k", linewidths=0.3, label=f"outlier (n={int(out.sum())})",
            rasterized=True, zorder=5,
        )
    lim = max(hg_log[mask].max(), kf_log[mask].max()) if mask.any() else 1
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.6)
    r = pearson(hg_log[mask], kf_log[mask]) if mask.sum() > 1 else float("nan")
    ax.set_xlabel("hg38 log2(TPM+1)")
    ax.set_ylabel("KOLF log2(TPM+1), raw summed")
    ax.set_title(f"{sample}\nPearson r = {r:.4f}")
    ax.legend(fontsize=7, markerscale=2)
    ax.set_aspect("equal", adjustable="box")

    ax2 = axes[1]
    mean_log_full = (hg_log + kf_log) / 2.0
    if inlier.any():
        ax2.scatter(mean_log_full[inlier], diff_log_all[inlier], s=1, alpha=0.15, c="#b2182b", rasterized=True)
    if out.any():
        ax2.scatter(
            mean_log_full[out], diff_log_all[out], s=28, alpha=0.85, c="#ff7f00",
            edgecolors="k", linewidths=0.3, rasterized=True, zorder=5,
        )
    ax2.axhline(0, color="k", lw=0.8, alpha=0.6)
    ax2.set_xlabel("mean log2(TPM+1)")
    ax2.set_ylabel("KOLF - hg38 log2(TPM+1)")
    ax2.set_title(f"Bland-Altman (bias={md:.3f}, outliers={int(outliers.sum())})")

    label_idx = np.where(outliers)[0]
    label_idx = label_idx[np.argsort(np.abs(diff_log_all[label_idx]))[::-1][:OUTLIER_LABEL_TOP]]
    for i in label_idx:
        name = pair_rows[i]["gene_name"]
        label = name if not name.startswith("ENSG") else pair_rows[i]["gene_id"]
        ax.annotate(label, (hg_log[i], kf_log[i]), fontsize=5, alpha=0.9, xytext=(3, 3), textcoords="offset points")
        ax2.annotate(label, (mean_log_full[i], diff_log_all[i]), fontsize=5, alpha=0.9, xytext=(3, 3), textcoords="offset points")

    fig.tight_layout()
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTDIR / f"{sample}_tpm_scatter_ba.png", dpi=150)
    plt.close(fig)


def analyze_pair(sample: str, hg: dict[str, dict], kf: dict[str, dict]) -> dict:
    genes = sorted(set(hg) & set(kf))
    hg_tpm = np.array([float(hg[g]["tpm"]) for g in genes])
    kf_raw = np.array([kolf_tpm_for_detection(kf, g) for g in genes])

    hg_log = np.log2(hg_tpm + 1.0)
    kf_log = np.log2(kf_raw + 1.0)
    expressed = (hg_tpm > 0) | (kf_raw > 0)
    for_outlier = expressed & (np.maximum(hg_tpm, kf_raw) >= OUTLIER_MIN_TPM)
    diff_log_all = kf_log - hg_log
    outliers = outlier_mask(diff_log_all, for_outlier)
    delta = kf_raw - hg_tpm

    pair_rows = [
        {
            "gene_id": g,
            "gene_name": hg[g]["gene_name"],
            "hg38_tpm": float(hg_tpm[i]),
            "kolf_tpm_raw": float(kf_raw[i]),
            "delta_tpm_raw": float(delta[i]),
            "log2fc_raw": float(diff_log_all[i]),
            "is_outlier": int(outliers[i]),
            "abs_delta_tpm": abs(float(delta[i])),
        }
        for i, g in enumerate(genes)
    ]
    _write_pair_tables(sample, pair_rows)

    mask = expressed
    diff_log = diff_log_all[mask]
    md = float(np.mean(diff_log)) if len(diff_log) else float("nan")
    _plot_pair_scatter(sample, hg_log, kf_log, diff_log_all, expressed, outliers, pair_rows, md)

    gained, lost = gain_loss(genes, hg, kf)
    row: dict = {
        "sample": sample,
        "n_genes": len(genes),
        "n_expressed": int(expressed.sum()),
        "pearson_log2tpm_expressed": pearson(hg_log[mask], kf_log[mask]) if mask.sum() > 1 else float("nan"),
        "spearman_log2tpm_expressed": spearman(hg_log[mask], kf_log[mask]) if mask.sum() > 1 else float("nan"),
        "mean_abs_delta_tpm_raw": float(np.mean(np.abs(delta))) if len(delta) else float("nan"),
        "bland_altman_bias_log2": md,
        "n_outliers": int(outliers.sum()),
        "genes_gained": gained,
        "genes_lost": lost,
    }
    for t in TPM_THRESHOLDS:
        key = f"{t}".replace(".", "p")
        row[f"hg38_genes_tpm_gt_{key}"] = count_genes_above(hg, t)
        row[f"kolf_genes_tpm_gt_{key}"] = count_genes_above(kf, t, kolf_raw=True)
        row[f"delta_genes_tpm_gt_{key}"] = row[f"kolf_genes_tpm_gt_{key}"] - row[f"hg38_genes_tpm_gt_{key}"]
    return row


def summary_figures(rows: list[dict]) -> None:
    if not rows:
        return
    samples = [r["sample"] for r in rows]
    x = np.arange(len(samples))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.bar(x, [r["pearson_log2tpm_expressed"] for r in rows], color="#4393c3")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Pearson r")
    ax.set_title("TPM concordance (raw summed KOLF)")
    ax.set_ylim(0, 1.05)

    ax = axes[0, 1]
    ax.bar(x - 0.2, [r["hg38_genes_tpm_gt_1p0"] for r in rows], 0.4, label="hg38", color="#92c5de")
    ax.bar(x + 0.2, [r["kolf_genes_tpm_gt_1p0"] for r in rows], 0.4, label="KOLF raw", color="#2166ac")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("# genes TPM > 1")
    ax.set_title("Gene detection (raw summed KOLF TPM)")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.bar(x - 0.2, [r["genes_gained"] for r in rows], 0.4, label="gained", color="#1b7837")
    ax.bar(x + 0.2, [r["genes_lost"] for r in rows], 0.4, label="lost", color="#d73027")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("# genes")
    ax.set_title("Gain/loss (TPM > 1 vs < 0.1)")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.bar(x, [r["n_outliers"] for r in rows], color="#ff7f00")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("# outliers")
    ax.set_title("TPM outliers (Bland-Altman)")

    fig.tight_layout()
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTDIR / "summary_across_samples.png", dpi=150)
    plt.close(fig)


def write_report(rows: list[dict], incomplete: list[str], excluded: list[str]) -> None:
    lines = [
        "# Reference comparison: KOLF diploid vs hg38",
        "",
        "StringTie reference-guided quantification (-e -L).",
        "",
        "**KOLF TPM:** sum across haplotypes (raw summed TPM for plots, correlation, outliers).",
        f"Outliers: Bland-Altman log2(TPM+1), TPM>={OUTLIER_MIN_TPM}, "
        f"|z|>{OUTLIER_MAD_Z} (MAD) and |log2FC|>{OUTLIER_MIN_ABS_LOG2FC}.",
        "",
        f"Complete pairs: {len(rows)}",
        f"Incomplete: {len(incomplete)}",
        f"Excluded samples: {len(excluded)}",
        "",
    ]
    if excluded:
        lines.extend(["Excluded (see excluded_samples.txt):"] + [f"- {s}" for s in excluded] + [""])
    if incomplete:
        lines.extend(["Incomplete samples:"] + [f"- {s}" for s in incomplete] + [""])

    if rows:
        mean_r = sum(r["pearson_log2tpm_expressed"] for r in rows) / len(rows)
        lines.extend(["## Summary", f"- Mean Pearson r (raw summed KOLF): {mean_r:.4f}", "", "## Per sample"])
        for r in rows:
            lines.append(
                f"- **{r['sample']}**: r={r['pearson_log2tpm_expressed']:.4f}, "
                f"TPM>1 hg38={r['hg38_genes_tpm_gt_1p0']} kolf_raw={r['kolf_genes_tpm_gt_1p0']} "
                f"(delta={r['delta_genes_tpm_gt_1p0']:+d}), "
                f"gain/loss={r['genes_gained']}/{r['genes_lost']}"
            )

    (OUTDIR / "REPORT.md").write_text("\n".join(lines) + "\n")


def write_outlier_summary(rows: list[dict]) -> None:
    if not rows:
        return
    pc_genes = load_named_protein_coding_genes()
    combined_path = OUTDIR / "tpm_outliers_all.tsv"
    pc_calls_path = OUTDIR / "tpm_outliers_protein_coding.tsv"
    strict_genes_path = OUTDIR / "tpm_outlier_genes_strict.tsv"

    all_calls: list[dict[str, str]] = []
    with combined_path.open("w") as out_fh:
        out_fh.write("sample\tgene_id\tgene_name\thg38_tpm\tkolf_tpm_raw\tdelta_tpm_raw\tlog2fc_raw\n")
        for row in rows:
            sample = row["sample"]
            path = PAIRDIR / f"{sample}_tpm_outliers.tsv"
            if not path.exists():
                continue
            with path.open() as in_fh:
                for gene in csv.DictReader(in_fh, delimiter="\t"):
                    rec = dict(gene)
                    rec["sample"] = sample
                    all_calls.append(rec)
                    out_fh.write(
                        f"{sample}\t{gene['gene_id']}\t{gene['gene_name']}\t"
                        f"{gene['hg38_tpm']}\t{gene['kolf_tpm_raw']}\t"
                        f"{gene['delta_tpm_raw']}\t{gene['log2fc_raw']}\n"
                    )

    pc_calls = [r for r in all_calls if is_named_protein_coding(r["gene_id"], r["gene_name"], pc_genes)]
    with pc_calls_path.open("w") as fh:
        fh.write("sample\tgene_id\tgene_name\thg38_tpm\tkolf_tpm_raw\tdelta_tpm_raw\tlog2fc_raw\n")
        for r in pc_calls:
            fh.write(
                f"{r['sample']}\t{r['gene_id']}\t{r['gene_name']}\t"
                f"{r['hg38_tpm']}\t{r['kolf_tpm_raw']}\t{r['delta_tpm_raw']}\t{r['log2fc_raw']}\n"
            )

    by_gene: dict[str, dict] = {}
    for r in pc_calls:
        gid = r["gene_id"]
        if gid not in by_gene:
            by_gene[gid] = {
                "gene_id": gid,
                "gene_name": pc_genes.get(gene_base_id(gid), r["gene_name"]),
                "samples": [],
                "log2fc_raw": [],
                "delta_tpm_raw": [],
            }
        g = by_gene[gid]
        g["samples"].append(r["sample"])
        g["log2fc_raw"].append(float(r["log2fc_raw"]))
        g["delta_tpm_raw"].append(float(r["delta_tpm_raw"]))

    agg = []
    for g in by_gene.values():
        n = len(g["samples"])
        mean_lfc = sum(g["log2fc_raw"]) / n
        agg.append(
            {
                "gene_id": g["gene_id"],
                "gene_name": g["gene_name"],
                "n_samples_outlier": n,
                "direction_consensus": "KOLF_higher" if mean_lfc > 0 else "KOLF_lower",
                "mean_log2fc_raw": mean_lfc,
                "mean_abs_log2fc_raw": sum(abs(x) for x in g["log2fc_raw"]) / n,
                "mean_delta_tpm_raw": sum(g["delta_tpm_raw"]) / n,
                "samples": ";".join(sorted(g["samples"])),
            }
        )
    agg.sort(key=lambda r: (-r["n_samples_outlier"], -r["mean_abs_log2fc_raw"]))

    cols = [
        "gene_id", "gene_name", "n_samples_outlier", "direction_consensus",
        "mean_log2fc_raw", "mean_abs_log2fc_raw", "mean_delta_tpm_raw", "samples",
    ]
    str_cols = {"gene_id", "gene_name", "direction_consensus", "samples", "n_samples_outlier"}
    with strict_genes_path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in agg:
            fh.write("\t".join(str(r[c]) if c in str_cols else f"{r[c]:.6f}" for c in cols) + "\n")

    print(
        f"Outlier genes: {len(all_calls)} calls, {len(pc_calls)} protein-coding calls, "
        f"{len(agg)} unique named protein-coding genes -> {strict_genes_path.name}"
    )


def run_tpm_comparison() -> None:
    exclude_patterns = load_excluded_patterns()
    kolf_id_map = load_kolf_id_map()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    by_sample: dict[str, dict[str, dict]] = {}
    all_samples: set[str] = set()
    excluded: list[str] = []

    for abund in sorted(QUANT.glob("*/gene_abundances.tab")):
        bam_base = abund.parent.name
        try:
            sample = sample_name(bam_base)
        except ValueError:
            continue  # unknown reference prefix (e.g. CHM13)
        if is_excluded(sample, exclude_patterns):
            if sample not in excluded:
                excluded.append(sample)
            continue
        all_samples.add(sample)
        if bam_base.startswith("hg38_"):
            by_sample.setdefault(sample, {})["hg38"] = load_hg38_gene_tpm(abund)
        elif bam_base.startswith("KOLF2.1Jv1.1_"):
            by_sample.setdefault(sample, {})["kolf"] = load_kolf_gene_tpm(abund, kolf_id_map)

    rows, incomplete = [], []
    for sample in sorted(all_samples):
        refs = by_sample.get(sample, {})
        if "hg38" not in refs or "kolf" not in refs:
            incomplete.append(sample)
            continue
        rows.append(analyze_pair(sample, refs["hg38"], refs["kolf"]))

    if rows:
        cols = list(rows[0].keys())
        with (OUTDIR / "summary.tsv").open("w") as fh:
            fh.write("\t".join(cols) + "\n")
            for row in rows:
                fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    summary_figures(rows)
    write_report(rows, incomplete, sorted(excluded))
    write_outlier_summary(rows)
    print(f"Analyzed {len(rows)} pairs (raw summed KOLF TPM) -> {OUTDIR}")
    if incomplete:
        print(f"Incomplete: {', '.join(incomplete)}")


def load_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short_sample(name: str) -> str:
    return name.replace("_cDNA_", "_").replace("_HiFi", "")


def load_pair_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hg, kf, is_outlier = [], [], []
    for row in csv.DictReader(path.open(), delimiter="\t"):
        h = float(row["hg38_tpm"])
        k = float(row["kolf_tpm_raw"])
        if h <= 0 and k <= 0:
            continue
        hg.append(h)
        kf.append(k)
        is_outlier.append(int(row.get("is_outlier", 0)) == 1)
    return np.array(hg, dtype=float), np.array(kf, dtype=float), np.array(is_outlier, dtype=bool)


def plot_scatter_across_samples(rows: list[dict]) -> None:
    n = len(rows)
    if n == 0:
        return
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), squeeze=False)

    global_lim = 0.0
    pair_data: list[tuple | None] = []
    for r in rows:
        path = PAIRDIR / f"{r['sample']}_tpm_pairs.tsv"
        if not path.exists():
            pair_data.append(None)
            continue
        hg, kf, is_outlier = load_pair_table(path)
        hg_log = np.log2(hg + 1.0)
        kf_log = np.log2(kf + 1.0)
        global_lim = max(global_lim, float(hg_log.max()), float(kf_log.max()))
        pair_data.append((hg_log, kf_log, is_outlier, r))

    lim = global_lim + (global_lim * 0.02 if global_lim else 0.1)

    for idx, r in enumerate(rows):
        ax = axes[idx // ncols][idx % ncols]
        data = pair_data[idx]
        if data is None:
            ax.set_visible(False)
            continue
        hg_log, kf_log, is_outlier, meta = data
        inlier = ~is_outlier
        if inlier.any():
            ax.scatter(hg_log[inlier], kf_log[inlier], s=0.8, alpha=0.15, c="#2166ac", rasterized=True)
        if is_outlier.any():
            ax.scatter(hg_log[is_outlier], kf_log[is_outlier], s=8, alpha=0.7, c="#ff7f00", rasterized=True)
        ax.plot([0, lim], [0, lim], "k--", lw=0.6, alpha=0.5)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{short_sample(meta['sample'])}\nr={meta['pearson_r']:.3f}", fontsize=8)
        if idx // ncols == nrows - 1:
            ax.set_xlabel("hg38 log2(TPM+1)", fontsize=7)
        else:
            ax.set_xticklabels([])
        if idx % ncols == 0:
            ax.set_ylabel("KOLF log2(TPM+1)", fontsize=7)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=6)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("TPM concordance across samples (raw summed KOLF vs hg38)", fontsize=11, y=1.01)
    fig.tight_layout()
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTDIR / "tpm_scatter_across_samples.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def summarize_tpm() -> None:
    raw = load_rows(SUMMARY_IN)
    if not raw:
        print(f"No TPM comparison summary yet: {SUMMARY_IN}")
        return

    rows = [
        {
            "sample": r["sample"],
            "n_genes": int(r["n_genes"]),
            "n_expressed": int(r["n_expressed"]),
            "pearson_r": float(r["pearson_log2tpm_expressed"]),
            "spearman_r": float(r["spearman_log2tpm_expressed"]),
            "bland_altman_bias_log2": float(r["bland_altman_bias_log2"]),
            "mean_abs_delta_tpm": float(r["mean_abs_delta_tpm_raw"]),
            "n_outliers": int(r["n_outliers"]),
            "hg38_genes_tpm_gt_1": int(r["hg38_genes_tpm_gt_1p0"]),
            "kolf_genes_tpm_gt_1": int(r["kolf_genes_tpm_gt_1p0"]),
            "delta_genes_tpm_gt_1": int(r["delta_genes_tpm_gt_1p0"]),
            "hg38_genes_tpm_gt_5": int(r["hg38_genes_tpm_gt_5p0"]),
            "kolf_genes_tpm_gt_5": int(r["kolf_genes_tpm_gt_5p0"]),
            "delta_genes_tpm_gt_5": int(r["delta_genes_tpm_gt_5p0"]),
            "genes_gained": int(r["genes_gained"]),
            "genes_lost": int(r["genes_lost"]),
        }
        for r in raw
    ]

    cols = list(rows[0].keys())
    with SUMMARY_OUT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            parts = []
            for c in cols:
                v = row[c]
                if isinstance(v, str):
                    parts.append(v)
                elif isinstance(v, int):
                    parts.append(str(v))
                else:
                    parts.append(f"{v:.6g}")
            fh.write("\t".join(parts) + "\n")

    samples = [short_sample(r["sample"]) for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    ax = axes[0, 0]
    ax.bar(x, [r["pearson_r"] for r in rows], color="#4393c3")
    ax.axhline(np.median([r["pearson_r"] for r in rows]), color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Pearson r")
    ax.set_title("TPM concordance (raw summed KOLF vs hg38)")
    ax.set_ylim(0.9, 1.0)

    ax = axes[0, 1]
    d_det = [r["delta_genes_tpm_gt_1"] for r in rows]
    colors = ["#1b7837" if d > 0 else "#d73027" if d < 0 else "#999999" for d in d_det]
    ax.bar(x, d_det, color=colors)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("delta genes detected (TPM > 1)")
    ax.set_title("Gene detection (KOLF raw - hg38)")

    ax = axes[1, 0]
    bias = [r["bland_altman_bias_log2"] for r in rows]
    colors = ["#1b7837" if b > 0 else "#d73027" for b in bias]
    ax.bar(x, bias, color=colors)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Bland-Altman bias (log2 TPM)")
    ax.set_title("Systematic TPM shift (KOLF - hg38)")

    ax = axes[1, 1]
    w = 0.35
    ax.bar(x - w / 2, [r["genes_gained"] for r in rows], w, label="gained", color="#1b7837")
    ax.bar(x + w / 2, [r["genes_lost"] for r in rows], w, label="lost", color="#d73027")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("# genes")
    ax.set_title("Gain/loss (TPM > 1 vs < 0.1)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOTDIR / "tpm_summary_across_samples.png", dpi=150)
    plt.close(fig)

    plot_scatter_across_samples(rows)

    med_r = np.median([r["pearson_r"] for r in rows])
    med_det = np.median(d_det)
    print(f"Wrote {SUMMARY_OUT}, tpm_summary_across_samples.png, tpm_scatter_across_samples.png ({len(rows)} samples)")
    print(f"  median Pearson r: {med_r:.4f}")
    print(f"  median delta genes TPM>1: {med_det:+.0f}")
