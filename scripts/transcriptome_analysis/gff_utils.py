from __future__ import annotations

import csv
import re
from pathlib import Path

from config import EXCLUDE_FILE, HG38_GFF, KOLF_ID_MAP, REF_PREFIXES


def load_excluded_patterns(path: Path = EXCLUDE_FILE) -> list[str]:
    if not path.exists():
        return []
    patterns: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        patterns.append(line)
    return patterns


def is_excluded(sample: str, patterns: list[str] | None = None) -> bool:
    pats = patterns if patterns is not None else load_excluded_patterns()
    return any(pat in sample for pat in pats)


def gene_base_id(gene_id: str) -> str:
    return re.sub(r"\.\d+(?:_\d+)?$", "", gene_id.strip())


def is_unnamed_gene(gene_name: str) -> bool:
    return bool(re.match(r"^ENSG", gene_name.strip()))


def sample_name(bam: str) -> str:
    name = bam[:-4] if bam.endswith(".bam") else bam
    for prefix in REF_PREFIXES.values():
        if name.startswith(prefix):
            return name[len(prefix):]
    raise ValueError(bam)


def ref_key(bam: str) -> str:
    for key, prefix in REF_PREFIXES.items():
        if bam.startswith(prefix):
            return key
    raise ValueError(bam)


def load_named_protein_coding_genes(gff_path: Path = HG38_GFF) -> dict[str, str]:
    out: dict[str, str] = {}
    if not gff_path.exists():
        return out
    for line in gff_path.open():
        if line.startswith("#") or "\tgene\t" not in line:
            continue
        attrs = line.rstrip().split("\t", 8)[-1]
        fields: dict[str, str] = {}
        for part in attrs.split(";"):
            if "=" not in part:
                continue
            key, val = part.split("=", 1)
            fields[key] = val
        if fields.get("gene_type") != "protein_coding":
            continue
        gid = gene_base_id(fields.get("gene_id", ""))
        name = fields.get("gene_name", "").strip()
        if not gid or is_unnamed_gene(name):
            continue
        out[gid] = name
    return out


def is_named_protein_coding(gene_id: str, gene_name: str, pc_genes: dict[str, str]) -> bool:
    base = gene_base_id(gene_id)
    if base not in pc_genes:
        return False
    if is_unnamed_gene(gene_name):
        return False
    return True


def load_hg38_gene_tpm(abund_path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with abund_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = gene_base_id(row.get("Gene ID") or "")
            if not gid:
                continue
            tpm = float(row.get("TPM") or 0.0)
            name = row.get("Gene Name") or gid
            if gid not in out:
                out[gid] = {"tpm": 0.0, "gene_name": name}
            out[gid]["tpm"] = float(out[gid]["tpm"]) + tpm
    return out


def load_kolf_id_map(path: Path = KOLF_ID_MAP) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cid = (row.get("cat_gene_id") or "").strip()
            if not cid:
                continue
            out[cid] = {
                "ensembl_base": (row.get("ensembl_base") or "").strip(),
                "gene_name": (row.get("gene_name") or "").strip(),
            }
    return out


def load_kolf_gene_tpm(abund_path: Path, id_map: dict[str, dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with abund_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw_id = (row.get("Gene ID") or "").strip()
            if not raw_id:
                continue
            meta = id_map.get(raw_id)
            if meta and meta["ensembl_base"]:
                gid = gene_base_id(meta["ensembl_base"])
                name = meta["gene_name"] or row.get("Gene Name") or gid
            else:
                gid = gene_base_id(raw_id)
                name = row.get("Gene Name") or gid
            if not gid:
                continue
            tpm = float(row.get("TPM") or 0.0)
            if gid not in out:
                out[gid] = {"tpm_raw": 0.0, "gene_name": name}
            out[gid]["tpm_raw"] = float(out[gid]["tpm_raw"]) + tpm
    return out
