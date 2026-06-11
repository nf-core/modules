#!/usr/bin/env python3
# Copyright (c) nf-core
# This software is licensed under the MIT License.
# SPDX-License-Identifier: MIT

"""Run rpbp's `get_orfs` reference-prep step.

Seeds the index layout (`chrName.txt`, the `transcript-index/` subdir) and
then calls the rpbp Python API, standing in for rpbp's `prepare-rpbp-genome`
umbrella script (which would also build bowtie2/STAR indices that are not
consumed here).
"""

import argparse
import gzip
import os
import platform
import shutil

import rpbp
import yaml
from pbiotools.misc import logging_utils
from rpbp.reference_preprocessing.prepare_rpbp_genome import get_orfs


def chr_names_from_fasta(fasta_path: str, out_path: str) -> None:
    """Write one chromosome name per line from a (possibly gzipped) FASTA."""
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith(">"):
                # FASTA header: drop leading '>' and trim at first whitespace.
                parts = line[1:].split()
                if parts:
                    out.write(parts[0] + "\\n")


prefix = "${prefix}"
fasta = "${fasta}"
gtf = "${gtf}"

base_path = "_rpbp_ref"
star_index = os.path.join(base_path, "star")
os.makedirs(os.path.join(base_path, "transcript-index"), exist_ok=True)
os.makedirs(star_index, exist_ok=True)
chr_names_from_fasta(fasta, os.path.join(star_index, "chrName.txt"))

config = {
    "genome_base_path": base_path,
    "genome_name": prefix,
    "fasta": fasta,
    "star_index": star_index,
}

# get_orfs reads rpbp's standard logging options off `args`; build them from
# rpbp's own parser, then set the few execution flags it also consults.
parser = argparse.ArgumentParser()
logging_utils.add_logging_options(parser)
rpbp_args = parser.parse_args([])
rpbp_args.do_not_call = False
rpbp_args.overwrite = False
rpbp_args.num_cpus = int("${task.cpus}")

get_orfs(gtf, rpbp_args, config, is_annotated=True, is_de_novo=False)

for src in [
    os.path.join(base_path, f"{prefix}.annotated.bed.gz"),
    os.path.join(base_path, "transcript-index", f"{prefix}.orfs-genomic.annotated.bed.gz"),
    os.path.join(base_path, "transcript-index", f"{prefix}.orfs-exons.annotated.bed.gz"),
]:
    shutil.move(src, os.path.basename(src))
shutil.rmtree(base_path)

with open("versions.yml", "w") as f:
    yaml.safe_dump(
        {"${task.process}": {"python": platform.python_version(), "rpbp": rpbp.__version__}},
        f,
    )
