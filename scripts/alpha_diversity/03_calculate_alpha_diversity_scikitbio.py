#!/usr/bin/env python3
"""Calculate alpha-diversity metrics with scikit-bio."""

from __future__ import annotations

import importlib.metadata
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from skbio.diversity.alpha import ace, chao1, observed_features, shannon, simpson


def is_count_like(values: np.ndarray) -> bool:
    finite = values[np.isfinite(values)]
    if finite.size == 0 or np.any(finite < 0):
        return False
    return bool(np.all(np.abs(finite - np.rint(finite)) < 1e-8) and np.nanmax(finite) > 1)


def main() -> None:
    output_dir = Path(os.environ.get("BLACKICE_ALPHA_OUTPUT_DIR", "results/alpha_diversity"))
    tables_dir = output_dir / "tables"
    logs_dir = output_dir / "logs"
    matrix_path = tables_dir / "alpha_diversity_abundance_matrix.csv"
    diagnostic_path = tables_dir / "alpha_diversity_diagnostic_summary.csv"

    if not matrix_path.exists():
        raise SystemExit("Prepared abundance matrix not found. Run 01_prepare_alpha_diversity_inputs.R first.")

    abundance = pd.read_csv(matrix_path)
    if "sampleID" not in abundance.columns:
        raise SystemExit("Prepared abundance matrix must contain sampleID.")

    sample_ids = abundance["sampleID"].astype(str)
    matrix = abundance.drop(columns=["sampleID"]).apply(pd.to_numeric, errors="raise")
    values = matrix.to_numpy(dtype=float)
    count_like = is_count_like(values)

    scale_decision = "unknown"
    if diagnostic_path.exists():
        diagnostic = pd.read_csv(diagnostic_path)
        hit = diagnostic.loc[diagnostic["item"] == "abundance_scale_decision", "value"]
        if not hit.empty:
            scale_decision = str(hit.iloc[0])

    rows = []
    for sample_id, row in zip(sample_ids, values):
        clean = np.nan_to_num(row, nan=0.0)
        if np.any(clean < 0):
            raise SystemExit(f"Negative abundance value found for sample {sample_id}.")
        if count_like:
            count_row = np.rint(clean).astype(int)
            chao_value = float(chao1(count_row))
            ace_value = float(ace(count_row))
            decision = "calculated"
        else:
            chao_value = math.nan
            ace_value = math.nan
            decision = "not_calculated_relative_or_non_count_input"

        rows.append(
            {
                "sampleID": sample_id,
                "Shannon": float(shannon(clean)),
                "Simpson": float(simpson(clean)),
                "Observed_richness": int(observed_features(clean)),
                "Chao1": chao_value,
                "ACE": ace_value,
                "tool": "scikit-bio",
                "abundance_scale_decision": scale_decision,
                "chao1_ace_decision": decision,
            }
        )

    result = pd.DataFrame(rows)
    result.to_csv(tables_dir / "alpha_diversity_scikitbio.csv", index=False)

    logs_dir.mkdir(parents=True, exist_ok=True)
    versions = {}
    for package in ("python", "pandas", "numpy", "scikit-bio"):
        if package == "python":
            import sys

            versions[package] = sys.version
        else:
            versions[package] = importlib.metadata.version(package)
    (logs_dir / "alpha_diversity_scikitbio_versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    print(f"Wrote scikit-bio alpha-diversity table: {tables_dir / 'alpha_diversity_scikitbio.csv'}")
    if not count_like:
        print("WARNING: Chao1 and ACE were not calculated because the abundance table is not count-like.")


if __name__ == "__main__":
    main()
