#!/usr/bin/env python3
# Written by Jonathan Manning (@pinin4fjords). Released under the MIT license.
"""Wrap rpbp's `get_periodic_lengths_and_offsets` helper.

Recreates the `rpbp_work/metagene-profiles/<prefix>-unique.periodic-offsets.csv.gz`
layout the helper expects, calls it with thresholds passed in from the
Nextflow process, and writes the resulting per-read-length offsets to a TSV.
"""

import os
import platform
import shutil

import pandas as pd
import rpbp
import yaml
from rpbp.ribo_utils.utils import get_periodic_lengths_and_offsets


def _parse_filter(value):
    """Parse a filter token to int/float or None."""
    if value == "None":
        return None
    try:
        return int(value)
    except ValueError:
        return float(value)


prefix = "${prefix}"
min_count = _parse_filter("${min_count}")
min_bf_mean = _parse_filter("${min_bf_mean}")
max_bf_var = _parse_filter("${max_bf_var}")
min_bf_lik = _parse_filter("${min_bf_lik}")

work_dir = os.path.join("rpbp_work", "metagene-profiles")
os.makedirs(work_dir, exist_ok=True)
shutil.copy(
    "${periodic_offsets}",
    os.path.join(work_dir, f"{prefix}-unique.periodic-offsets.csv.gz"),
)

config = dict(
    riboseq_data="rpbp_work",
    min_metagene_profile_count=min_count,
    min_metagene_bf_mean=min_bf_mean,
    max_metagene_bf_var=max_bf_var,
    min_metagene_bf_likelihood=min_bf_lik,
)

lengths, offsets = get_periodic_lengths_and_offsets(config, prefix, is_unique=True)
if len(lengths) == 0:
    raise SystemExit(
        "No periodic read lengths passed filters; "
        "check min_count/min_bf_mean thresholds and metagene Bayes-factor output."
    )
pd.DataFrame({"length": lengths, "offset": offsets}).to_csv(
    f"{prefix}.periodic_lengths_offsets.tsv",
    sep="\\t",
    index=False,
)

with open("versions.yml", "w") as f:
    yaml.safe_dump(
        {"${task.process}": {"python": platform.python_version(), "rpbp": rpbp.__version__}},
        f,
    )
