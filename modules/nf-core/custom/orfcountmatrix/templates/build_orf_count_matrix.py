#!/usr/bin/env python3
"""Pivot per-sample ORF P-site count TSVs into a single ORF x sample matrix.

Per-sample input is `sample_id<TAB>orf_id<TAB>count` (sample_id is identical
on every row, prepended upstream by the per-sample counter). Output rows
follow the BED12 catalogue's 4th column in catalogue order; ORFs absent
from a sample are zero-filled, and ORFs absent from every sample still
appear as a row of zeros so the matrix is keyed on the catalogue.
"""

import platform
from pathlib import Path

import pandas as pd
import yaml

counts = pd.concat(
    pd.read_csv(p, sep="\\t", names=["sample", "orf_id", "count"]) for p in sorted(Path("counts").iterdir())
)

catalogue_orfs = (
    pd.read_csv("$orf_catalogue_bed12", sep="\\t", header=None, comment="#", usecols=[3])[3].drop_duplicates().tolist()
)

matrix = (
    counts.pivot_table(index="orf_id", columns="sample", values="count", aggfunc="sum", fill_value=0)
    .reindex(catalogue_orfs, fill_value=0)
    .astype(int)
)
matrix.to_csv("${prefix}.tsv", sep="\\t")

with open("versions.yml", "w") as f:
    yaml.safe_dump(
        {"${task.process}": {"python": platform.python_version(), "pandas": pd.__version__}},
        f,
    )
