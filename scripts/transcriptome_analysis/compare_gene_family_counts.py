#!/usr/bin/env python3
"""Gene-family copy-number (CNV) comparison across annotations.
Requirements: python3 (pysam, biopython, networkx, numpy), blastp/makeblastdb.

Example
-------
    python compare_gene_family_counts.py \\
        --grch38-gff3   GRCh38.gff3 \\
        --grch38-fasta  GRCh38.fa \\
        --chm13-gff3    CHM13.gff3 \\
        --kolf-pat-gff3 KOLF_PAT_noY.gff3 \\
        --kolf-mat-gff3 KOLF_MAT.gff3 \\
        --out-dir       gene_family_cnv \\
        --threads       32
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pysam
from Bio.Seq import Seq

AUTOSOMES: frozenset[str] = frozenset(f"chr{i}" for i in range(1, 23))

ANNOTATIONS = ("CHM13", "GRCh38", "KOLF_PAT", "KOLF_MAT")
HAPLOTYPES = ("KOLF_PAT", "KOLF_MAT")


def parse_attributes(attr_field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def ensembl_base(gene_id: str) -> str:
    return re.sub(r"\.\d+(?:_\d+)?$", "", gene_id.strip())

def extract_mane_proteins(
    gff3_path: Path,
    fasta_path: Path,
    out_fasta: Path,
    *,
    autosomes_only: bool = True,
) -> dict[str, str]:
    allowed_chroms = AUTOSOMES if autosomes_only else None

    # First pass: pick the MANE Select transcript per protein-coding gene.
    mane_tx: dict[str, dict] = {}   # transcript_id -> {gene, name, chrom, strand}
    cds_by_tx: dict[str, list[tuple[int, int, int]]] = defaultdict(list)  # tx -> [(start,end,phase)]
    with gff3_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, chrom = fields[2], fields[0]
            if feature not in ("transcript", "mRNA", "CDS"):
                continue
            if allowed_chroms is not None and chrom not in allowed_chroms:
                continue
            attrs = parse_attributes(fields[8])
            if feature in ("transcript", "mRNA"):
                if "MANE_Select" not in attrs.get("tag", ""):
                    continue
                if attrs.get("gene_type") != "protein_coding":
                    continue
                gene = ensembl_base(attrs.get("gene_id", ""))
                tx = attrs.get("transcript_id", "")
                if not gene or not tx or gene in {v["gene"] for v in mane_tx.values()}:
                    continue
                mane_tx[tx] = {
                    "gene": gene,
                    "name": attrs.get("gene_name", gene),
                    "chrom": chrom,
                    "strand": fields[6],
                }
            else:  # CDS
                parent = attrs.get("Parent", "")
                if not parent:
                    continue
                phase = int(fields[7]) if fields[7].isdigit() else 0
                cds_by_tx[parent].append((int(fields[3]), int(fields[4]), phase))

    gene_symbol: dict[str, str] = {}
    fasta = pysam.FastaFile(str(fasta_path))
    n_written = 0
    with out_fasta.open("w") as out:
        for tx, meta in mane_tx.items():
            cds = sorted(cds_by_tx.get(tx, []))
            if not cds:
                continue
            seq = "".join(fasta.fetch(meta["chrom"], s - 1, e) for s, e, _ in cds)
            nt = Seq(seq)
            if meta["strand"] == "-":
                nt = nt.reverse_complement()
                first_phase = cds[-1][2]
            else:
                first_phase = cds[0][2]
            nt = nt[first_phase:]
            protein = str(nt.translate(to_stop=True))
            if len(protein) < 10:
                continue
            gene = meta["gene"]
            gene_symbol[gene] = meta["name"]
            out.write(f">{gene} {meta['name']}\n{protein}\n")
            n_written += 1
    fasta.close()
    print(f"[mane] wrote {n_written} MANE proteins -> {out_fasta.name}")
    return gene_symbol

def run_blastp_all_vs_all(
    protein_fasta: Path,
    blast_out: Path,
    *,
    threads: int,
    evalue: float,
    force: bool,
) -> None:
    if blast_out.exists() and not force:
        print(f"[blast] skip (exists): {blast_out.name}")
        return
    db_prefix = protein_fasta.with_suffix(".blastdb")
    print("[blast] makeblastdb")
    subprocess.run(
        ["makeblastdb", "-in", str(protein_fasta), "-dbtype", "prot",
         "-out", str(db_prefix)],
        check=True, stdout=subprocess.DEVNULL,
    )
    print(f"[blast] blastp all-vs-all ({threads} threads)")
    subprocess.run(
        ["blastp", "-query", str(protein_fasta), "-db", str(db_prefix),
         "-out", str(blast_out), "-evalue", str(evalue),
         "-num_threads", str(threads), "-max_target_seqs", "500",
         "-outfmt", "6 qseqid sseqid pident length qlen slen"],
        check=True,
    )
    print(f"[blast] wrote {blast_out.name}")


def cluster_proteins(
    blast_out: Path,
    genes: set[str],
    *,
    min_identity: float,
    min_coverage: float,
) -> dict[str, int]:
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    with blast_out.open() as handle:
        for line in handle:
            q, s, pident, length, qlen, slen = line.rstrip("\n").split("\t")
            if q == s or q not in genes or s not in genes:
                continue
            if float(pident) < min_identity:
                continue
            if int(length) / min(int(qlen), int(slen)) < min_coverage:
                continue
            graph.add_edge(q, s)

    gene_cluster: dict[str, int] = {}
    for cid, component in enumerate(nx.connected_components(graph)):
        for gene in component:
            gene_cluster[gene] = cid
    n_multi = sum(1 for c in nx.connected_components(graph) if len(c) > 1)
    print(f"[cluster] {len(set(gene_cluster.values()))} families "
          f"({n_multi} multi-gene) from {len(genes)} genes")
    return gene_cluster

def count_family_copies(
    gff3_path: Path,
    gene_cluster: dict[str, int],
    *,
    is_grch38: bool,
) -> dict[int, int]:
    source_key = "gene_id" if is_grch38 else "source_gene"
    counts: dict[int, int] = defaultdict(int)
    with gff3_path.open() as handle:
        for line in handle:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            if attrs.get("gene_biotype", attrs.get("gene_type")) != "protein_coding":
                continue
            source = ensembl_base(attrs.get(source_key, ""))
            cluster = gene_cluster.get(source)
            if cluster is None:
                continue
            counts[cluster] += 1
    return counts

def call_cnv_events(
    copies: dict[str, dict[int, int]],
    gene_cluster: dict[str, int],
    gene_symbol: dict[str, str],
) -> tuple[list[dict], dict[str, dict[str, float]]]:
    cluster_members: dict[int, list[str]] = defaultdict(list)
    for gene, cid in gene_cluster.items():
        cluster_members[cid].append(gene)

    totals = {
        hap: {"duplications": 0.0, "deletions": 0.0, "families": 0}
        for hap in HAPLOTYPES
    }
    rows: list[dict] = []
    for cid, members in cluster_members.items():
        vals = {ann: copies[ann].get(cid, 0) for ann in ANNOTATIONS}
        median = float(vals["CHM13"])
        row: dict = {
            "family_id": cid,
            "n_reference_genes": len(members),
            "example_symbols": ",".join(
                sorted({gene_symbol.get(g, g) for g in members})[:5]
            ),
            "baseline_copy_number": median,
            **{f"copies_{ann}": vals[ann] for ann in ANNOTATIONS},
        }
        event = False
        for hap in HAPLOTYPES:
            dev = vals[hap] - median
            row[f"{hap}_delta"] = dev
            if dev > 0:
                row[f"{hap}_event"] = "duplication"
                totals[hap]["duplications"] += dev
                totals[hap]["families"] += 1
                event = True
            elif dev < 0:
                row[f"{hap}_event"] = "deletion"
                totals[hap]["deletions"] += -dev
                totals[hap]["families"] += 1
                event = True
            else:
                row[f"{hap}_event"] = ""
        if event:
            rows.append(row)

    rows.sort(key=lambda r: (
        -max(abs(r[f"{h}_delta"]) for h in HAPLOTYPES),
        r["family_id"],
    ))
    return rows, totals


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--grch38-gff3", required=True, type=Path,
                        help="GRCh38 annotation with MANE_Select tags and CDS features.")
    parser.add_argument("--grch38-fasta", required=True, type=Path,
                        help="GRCh38 genome FASTA (indexed; .fai alongside).")
    parser.add_argument("--chm13-gff3", required=True, type=Path,
                        help="CHM13 CAT annotation GFF3.")
    parser.add_argument("--kolf-pat-gff3", required=True, type=Path,
                        help="KOLF paternal-haplotype CAT GFF3.")
    parser.add_argument("--kolf-mat-gff3", required=True, type=Path,
                        help="KOLF maternal-haplotype CAT GFF3.")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output/cache directory (created if missing).")
    parser.add_argument("--threads", type=int, default=8,
                        help="BLAST threads (default: 8).")
    parser.add_argument("--min-identity", type=float, default=90.0,
                        help="Min percent identity for a family edge (default: 90).")
    parser.add_argument("--min-coverage", type=float, default=0.8,
                        help="Min aligned fraction of the shorter protein (default: 0.8).")
    parser.add_argument("--evalue", type=float, default=1e-6,
                        help="BLAST e-value cutoff (default: 1e-6).")
    parser.add_argument("--include-sex-chromosomes", action="store_true",
                        help="Include chrX/chrY MANE genes (default: autosomes only).")
    parser.add_argument("--force", action="store_true",
                        help="Recompute cached intermediates (proteins, BLAST).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    protein_fasta = args.out_dir / "mane_proteins.faa"
    blast_out = args.out_dir / "mane_blastp.tsv"

    if protein_fasta.exists() and not args.force:
        print(f"[mane] skip (exists): {protein_fasta.name}")
        gene_symbol = {}
        for line in protein_fasta.read_text().splitlines():
            if not line.startswith(">"):
                continue
            parts = line[1:].split(None, 1)
            gene_symbol[parts[0]] = parts[1].strip() if len(parts) > 1 else parts[0]
    else:
        gene_symbol = extract_mane_proteins(
            args.grch38_gff3, args.grch38_fasta, protein_fasta,
            autosomes_only=not args.include_sex_chromosomes,
        )

    genes = set(gene_symbol)
    run_blastp_all_vs_all(
        protein_fasta, blast_out,
        threads=args.threads, evalue=args.evalue, force=args.force,
    )
    gene_cluster = cluster_proteins(
        blast_out, genes,
        min_identity=args.min_identity, min_coverage=args.min_coverage,
    )

    copies = {
        "GRCh38": count_family_copies(args.grch38_gff3, gene_cluster, is_grch38=True),
        "CHM13": count_family_copies(args.chm13_gff3, gene_cluster, is_grch38=False),
        "KOLF_PAT": count_family_copies(args.kolf_pat_gff3, gene_cluster, is_grch38=False),
        "KOLF_MAT": count_family_copies(args.kolf_mat_gff3, gene_cluster, is_grch38=False),
    }

    rows, totals = call_cnv_events(copies, gene_cluster, gene_symbol)

    detail_cols = (
        ["family_id", "n_reference_genes", "example_symbols", "baseline_copy_number"]
        + [f"copies_{ann}" for ann in ANNOTATIONS]
        + [c for hap in HAPLOTYPES for c in (f"{hap}_delta", f"{hap}_event")]
    )
    detail_path = args.out_dir / "gene_family_cnv_events.tsv"
    write_tsv(detail_path, detail_cols, rows)

    summary_rows = [
        {
            "haplotype": hap,
            "duplication_events": round(totals[hap]["duplications"]),
            "deletion_events": round(totals[hap]["deletions"]),
            "total_events": round(totals[hap]["duplications"] + totals[hap]["deletions"]),
            "families_affected": totals[hap]["families"],
        }
        for hap in HAPLOTYPES
    ]
    summary_path = args.out_dir / "gene_family_cnv_summary.tsv"
    write_tsv(
        summary_path,
        ["haplotype", "duplication_events", "deletion_events", "total_events",
         "families_affected"],
        summary_rows,
    )

    for row in summary_rows:
        print(f"{row['haplotype']}: {row['total_events']} CNV events "
              f"({row['duplication_events']} dup, {row['deletion_events']} del) "
              f"across {row['families_affected']} families")
    print(f"Details -> {detail_path}")
    print(f"Summary -> {summary_path}")


if __name__ == "__main__":
    main()
