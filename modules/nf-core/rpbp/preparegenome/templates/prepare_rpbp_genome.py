#!/usr/bin/env python3
"""Run rpbp's `get_orfs` reference-prep step.

Drops in for the original bash + heredoc combo that prepares the index
layout (`chrName.txt` for the STAR stub, the `transcript-index/`
subdir) and then calls the rpbp Python API.
"""

import argparse
import os
import platform
import shutil

import yaml

import rpbp
from rpbp.reference_preprocessing.prepare_rpbp_genome import get_orfs


def chr_names_from_fasta(fasta_path: str, out_path: str) -> None:
    """Write one chromosome name per line from a (possibly gzipped) FASTA."""
    import gzip

    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith(">"):
                # FASTA header: drop leading '>' and trim at first whitespace.
                out.write(line[1:].split()[0] + "\\n")


prefix = "${prefix}"
name   = "${name}"
fasta  = "${fasta}"
gtf    = "${gtf}"
ncpus  = int("${task.cpus}")

os.makedirs(os.path.join(prefix, "transcript-index"), exist_ok=True)
os.makedirs(os.path.join(prefix, "star"), exist_ok=True)

chr_names_from_fasta(fasta, os.path.join(prefix, "star", "chrName.txt"))

config = {
    "genome_base_path": prefix,
    "genome_name":      name,
    "gtf":              gtf,
    "fasta":            fasta,
    "star_index":       os.path.join(prefix, "star"),
}

args = argparse.Namespace(
    do_not_call=False,
    overwrite=False,
    num_cpus=ncpus,
    log_file="",
    enable_ext_logging=False,
    log_stdout=False,
    no_log_stderr=False,
    logging_level="WARNING",
    file_logging_level="NOTSET",
    stdout_logging_level="NOTSET",
    stderr_logging_level="NOTSET",
)

get_orfs(config["gtf"], args, config, is_annotated=True, is_de_novo=False)

shutil.rmtree(os.path.join(prefix, "star"), ignore_errors=True)

with open("versions.yml", "w") as f:
    yaml.safe_dump(
        {"${task.process}": {"python": platform.python_version(), "rpbp": rpbp.__version__}},
        f,
    )
