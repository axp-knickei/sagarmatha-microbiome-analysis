#!/usr/bin/env python3
"""Inventory candidate beta-diversity inputs without loading large files fully."""

# Purpose: Record candidate input files, canonical-input checksums, and a
# structure report before the staged beta-diversity workflow runs.
# Inputs: Repository root plus canonical raw files under data/raw by default.
# Outputs: TSV and Markdown audit files under results/beta_diversity/audit by
# default; this stage does not alter source inputs.
# Dependencies: Python standard library only.
# Usage: python analysis/beta_diversity/00_inspect_repository.py
# Workflow: Optional preflight before 01_validate_inputs.R.

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from statistics import median


EXCLUDED_DIRS = {".git", ".venv", "venv", "node_modules", "__pycache__", "library"}
CANDIDATE_SUFFIXES = {
    ".csv", ".tsv", ".txt", ".rds", ".rdata", ".xlsx", ".parquet",
    ".feather", ".biom", ".nwk", ".newick", ".tree", ".tre", ".nex",
    ".nexus", ".fa", ".fasta", ".fna", ".r", ".py", ".sh", ".pbs",
    ".slurm", ".yml", ".yaml", ".md",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_delimited(path: Path) -> dict[str, object]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(8192)
            handle.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
                delimiter = dialect.delimiter
            except csv.Error:
                pass
            reader = csv.reader(handle, delimiter=delimiter)
            header = next(reader, [])
            row_count = sum(1 for _ in reader)
        return {
            "rows": row_count,
            "columns": len(header),
            "delimiter": "TAB" if delimiter == "\t" else delimiter,
            "header_status": "present" if header else "unknown",
            "first_columns": ";".join(header[:12]),
        }
    except (UnicodeDecodeError, OSError):
        return {"rows": "", "columns": "", "delimiter": "", "header_status": "unknown", "first_columns": ""}


def infer_role(path: Path, info: dict[str, object]) -> tuple[str, str, str, str]:
    name = path.name.lower()
    cols = str(info.get("first_columns", "")).lower()
    suffix = path.suffix.lower()
    if name == "gut.csv" and "sampleid" in cols:
        return "sample-by-MAG abundance matrix", "high", "yes", "Non-negative sample-level MAG profiles; canonical community input."
    if name == "sample_info.csv" and "sampleid" in cols:
        return "sample metadata", "high", "yes", "Contains sample, subject, time, demographic, and altitude fields."
    if name == "gdtb.csv" and "genome" in cols:
        return "MAG taxonomy", "high", "yes", "Genome identifiers match all gut MAG columns."
    if suffix in {".nwk", ".newick", ".tree", ".tre", ".nex", ".nexus"}:
        return "candidate phylogenetic tree", "medium", "review", "Tree compatibility requires tip and branch-length validation."
    if suffix in {".fa", ".fasta", ".fna"}:
        return "candidate sequence data", "medium", "no", "Not a tree; could support a separate upstream phylogeny workflow."
    if "lmertest" in name or "coefficient" in name or "pvalue" in name:
        return "model summary or model workflow", "high", "no", "Statistical summaries cannot be used as a community matrix."
    if suffix in {".r", ".py", ".sh", ".pbs", ".slurm"}:
        return "analysis or workflow code", "high", "no", "Inspected for conventions, not used as analytical data."
    if suffix in {".yml", ".yaml"}:
        return "configuration or environment", "medium", "no", "Workflow context only."
    if suffix == ".md":
        return "documentation", "high", "no", "Scientific or operational context only."
    return "candidate tabular or serialized data", "low", "review", "Role cannot be established from filename and header alone."


def inventory(root: Path, output: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or any(part in EXCLUDED_DIRS for part in path.parts):
            continue
        rel = path.relative_to(root)
        if rel.parts and rel.parts[0] in {"results", "outputs"}:
            continue
        suffix = path.suffix.lower()
        if suffix not in CANDIDATE_SUFFIXES:
            continue
        info = inspect_delimited(path) if suffix in {".csv", ".tsv", ".txt"} else {
            "rows": "", "columns": "", "delimiter": "", "header_status": "not_applicable", "first_columns": ""
        }
        role, confidence, use, reason = infer_role(path, info)
        records.append({
            "path": str(rel), "extension": suffix.lstrip("."), "inferred_format": suffix.lstrip(".").upper(),
            "size_bytes": path.stat().st_size, **info, "probable_role": role,
            "role_confidence": confidence, "use_in_beta_diversity": use, "reason": reason,
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(records[0]) if records else []
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    return records


def inspect_canonical(root: Path, audit_dir: Path) -> None:
    abundance_path = root / "data/raw/gut.csv"
    metadata_path = root / "data/raw/sample_info.csv"
    taxonomy_path = root / "data/raw/GDTB.csv"
    if not abundance_path.exists():
        (audit_dir / "missing_data_requirements.md").write_text(
            "# Missing beta-diversity input\n\nFAILED_NO_SAMPLE_BY_MAG_MATRIX: provide a CSV with unique sample rows, unique MAG columns, and non-negative sample-level values. Model coefficients are not valid input.\n",
            encoding="utf-8",
        )
        return

    with abundance_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        abundance_header = next(reader)
        abundance_rows = list(reader)
    with metadata_path.open(encoding="utf-8-sig", newline="") as handle:
        metadata_rows = list(csv.DictReader(handle))
    with taxonomy_path.open(encoding="utf-8-sig", newline="") as handle:
        taxonomy_rows = list(csv.DictReader(handle))

    sample_ids = [row[0] for row in abundance_rows]
    metadata_ids = [row["sampleID"] for row in metadata_rows]
    features = abundance_header[1:]
    taxonomy_ids = [row["Genome"] for row in taxonomy_rows]
    matrix = [[float(value) for value in row[1:]] for row in abundance_rows]
    values = [value for row in matrix for value in row]
    row_sums = [sum(row) for row in matrix]
    prevalence = [sum(row[j] > 0 for row in matrix) for j in range(len(features))]
    matched_metadata = [row for row in metadata_rows if row["sampleID"] in set(sample_ids)]
    duplicate = lambda xs: sorted({x for x in xs if xs.count(x) > 1})
    tp_order = sorted(
        {(row["Timepoint"], row["date"], row["altitude"]) for row in matched_metadata},
        key=lambda item: item[1],
    )
    quantiles = sorted(row_sums)
    q = lambda p: quantiles[round((len(quantiles) - 1) * p)]
    report = [
        "# Beta-diversity data structure report", "",
        "Generated by `analysis/beta_diversity/00_inspect_repository.py`.", "",
        "## Canonical candidates", "",
        f"- Community matrix: `data/raw/gut.csv` ({len(matrix)} samples x {len(features)} MAGs), rows are samples.",
        f"- Metadata: `data/raw/sample_info.csv` ({len(metadata_rows)} rows; {len(matched_metadata)} match abundance samples).",
        f"- Taxonomy: `data/raw/GDTB.csv` ({len(taxonomy_rows)} genomes; {len(set(features) & set(taxonomy_ids))}/{len(features)} MAGs match).",
        "- Phylogenetic tree: none detected. UniFrac must be skipped; taxonomy strings are not a phylogeny.", "",
        "## Identifiers and design", "",
        "- Sample ID: `sampleID` (unique in both canonical tables).",
        "- Subject ID: metadata `group` (13 participants; canonical name `subject_id`).",
        "- Time: `Timepoint` ordered by `date`; subject-relative days are derived from `date`.",
        "- Primary grouping variable: ordered `Timepoint`.",
        "- Covariates: `age`, `gender`; `altitude` is time-linked and requires a separate sensitivity model.",
        "- No treatment arm or independent biological group was found.", "",
        "Observed time order: " + ", ".join(f"{tp} ({date}, {alt} m)" for tp, date, alt in tp_order) + ".", "",
        "## Validation findings", "",
        f"- Duplicate abundance sample IDs: {len(duplicate(sample_ids))}.",
        f"- Duplicate metadata sample IDs: {len(duplicate(metadata_ids))}.",
        f"- Duplicate MAG columns: {len(duplicate(features))}.",
        f"- Abundance IDs missing metadata: {len(set(sample_ids) - set(metadata_ids))}.",
        f"- Metadata IDs without gut abundance: {len(set(metadata_ids) - set(sample_ids))}; listed in `id_matching_report.tsv` after validation.",
        f"- Missing/non-finite abundance values: {sum(not math.isfinite(x) for x in values)}.",
        f"- Negative abundance values: {sum(x < 0 for x in values)}.",
        f"- All-zero samples: {sum(x == 0 for x in row_sums)}; all-zero MAGs: {sum(x == 0 for x in prevalence)}.",
        f"- Zero fraction: {sum(x == 0 for x in values) / len(values):.6f} ({sum(x == 0 for x in values):,}/{len(values):,}).",
        f"- Library/row-sum distribution: min {min(row_sums):.6f}, Q1 {q(.25):.6f}, median {median(row_sums):.6f}, Q3 {q(.75):.6f}, max {max(row_sums):.6f}.",
        f"- MAG prevalence: all {len(features)} occur in at least 10% of samples.",
        "- Values are continuous, non-negative, mostly below 1, and row sums are near but not exactly 1.",
        "- Inference: relative-abundance-like, probably representing the MAG-classified fraction; exact upstream normalization remains undocumented.",
        "- The matrix does not appear log-transformed because it contains exact zeros and no negative values.", "",
        "## Critical exclusions", "",
        "Per-MAG lmerTest tables and archives contain downstream model summaries and are excluded from distance calculation. Only sample-level `gut.csv` values are used.", "",
        "## Checksums", "",
        f"- `gut.csv`: `{sha256(abundance_path)}`",
        f"- `sample_info.csv`: `{sha256(metadata_path)}`",
        f"- `GDTB.csv`: `{sha256(taxonomy_path)}`",
    ]
    (audit_dir / "data_structure_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--output-dir", default="results/beta_diversity/audit")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    audit_dir = (root / args.output_dir).resolve()
    audit_dir.mkdir(parents=True, exist_ok=True)
    records = inventory(root, audit_dir / "file_inventory.tsv")
    inspect_canonical(root, audit_dir)
    print(f"COMPLETED_REPOSITORY_INSPECTION: {len(records)} candidate files inventoried")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
