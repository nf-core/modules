#!/usr/bin/env python3
# Written by Jonathan Manning (@pinin4fjords). Released under the MIT license.
"""Wrap rpbp's `get_periodic_lengths_and_offsets` helper.

Recreates the `rpbp_work/metagene-profiles/<prefix>-unique.periodic-offsets.csv.gz`
layout the helper expects, calls it with thresholds passed in from the
Nextflow process, and writes the resulting per-read-length offsets to a TSV.
"""

import argparse
import os
import shlex
import shutil

import pandas as pd
from rpbp.ribo_utils.utils import get_periodic_lengths_and_offsets


prefix = "${prefix}"

parser = argparse.ArgumentParser()
parser.add_argument("--min-count", type=int, default=1000, dest="min_count")
parser.add_argument("--min-bf-mean", type=float, default=5.0, dest="min_bf_mean")
parser.add_argument("--max-bf-var", type=float, default=None, dest="max_bf_var")
parser.add_argument("--min-bf-likelihood", type=float, default=0.5, dest="min_bf_likelihood")
args = parser.parse_args(shlex.split("${task_ext_args}"))

work_dir = os.path.join("rpbp_work", "metagene-profiles")
os.makedirs(work_dir, exist_ok=True)
shutil.copy(
    "${periodic_offsets}",
    os.path.join(work_dir, f"{prefix}-unique.periodic-offsets.csv.gz"),
)

config = dict(
    riboseq_data="rpbp_work",
    min_metagene_profile_count=args.min_count,
    min_metagene_bf_mean=args.min_bf_mean,
    max_metagene_bf_var=args.max_bf_var,
    min_metagene_bf_likelihood=args.min_bf_likelihood,
)

lengths, offsets = get_periodic_lengths_and_offsets(config, prefix, is_unique=True)
if len(lengths) == 0:
    raise SystemExit(
        "No periodic read lengths passed filters; "
        "check --min-count/--min-bf-mean thresholds and metagene Bayes-factor output."
    )
pd.DataFrame({"length": lengths, "offset": offsets}).to_csv(
    f"{prefix}.tsv",
    sep="\\t",
    index=False,
)
