#!/usr/bin/env python

# Combine featureCounts biotype counts with a header for MultiQC,
# then calculate feature percentages for the general stats table.
#
# Written by Senthilkumar Panneerselvam and released under the MIT license.
# Adapted for nf-core/modules by Jonathan Manning.

import logging
import os
import platform

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

# Template variables from Nextflow
count_file = "${count}"
header_file = "${header}"
prefix = "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${meta.id}"
sample_name = "${meta.id}"


def prepare_biotype_counts(count_file, header_file, prefix):
    """Cut gene ID and count columns from featureCounts output, prepend header."""
    outfile = f"{prefix}.biotype_counts_mqc.tsv"

    with open(header_file) as hf:
        header = hf.read()

    with open(count_file) as cf:
        lines = cf.readlines()[2:]  # skip featureCounts header lines

    with open(outfile, "w") as of:
        of.write(header)
        for line in lines:
            fields = line.strip().split("\\t")
            if len(fields) >= 7:
                of.write(f"{fields[0]}\\t{fields[6]}\\n")

    return outfile


mqc_main = """#id: 'biotype-gs'
#plot_type: 'generalstats'
#pconfig:"""

mqc_pconf = """#    percent_{ft}:
#        title: '% {ft}'
#        namespace: 'Biotype Counts'
#        description: '% reads overlapping {ft} features'
#        max: 100
#        min: 0
#        scale: 'RdYlGn-rev'"""


def mqc_feature_stat(bfile, features, outfile, sname=None):
    """Calculate features percentage for biotype counts."""
    if not sname:
        sname = os.path.splitext(os.path.basename(bfile))[0]

    fcounts = {}
    try:
        with open(bfile) as bfl:
            for ln in bfl:
                if ln.startswith("#"):
                    continue
                parts = ln.strip().split("\\t")
                if len(parts) == 2:
                    ft, cn = parts
                    fcounts[ft] = float(cn)
    except Exception:
        logger.error(f"Trouble reading the biocount file {bfile}")
        return

    total_count = sum(fcounts.values())
    if total_count == 0:
        logger.error("No biocounts found, exiting")
        return

    fpercent = {f: (fcounts[f] / total_count) * 100 if f in fcounts else 0 for f in features}
    if len(fpercent) == 0:
        logger.error(f"None of given features '{', '.join(features)}' found in the biocount file {bfile}")
        return

    out_head, out_value, out_mqc = ("Sample", f"'{sname}'", mqc_main)
    for ft, pt in fpercent.items():
        out_head = f"{out_head}\\tpercent_{ft}"
        out_value = f"{out_value}\\t{pt}"
        out_mqc = f"{out_mqc}\\n{mqc_pconf.format(ft=ft)}"

    with open(outfile, "w") as ofl:
        out_final = "\\n".join([out_mqc, out_head, out_value]).strip()
        ofl.write(out_final + "\\n")


if __name__ == "__main__":
    # Step 1: Prepare biotype counts TSV from featureCounts output
    biotype_file = prepare_biotype_counts(count_file, header_file, prefix)

    # Step 2: Calculate rRNA percentage for MultiQC general stats
    mqc_feature_stat(
        biotype_file,
        ["rRNA"],
        f"{prefix}.biotype_counts_rrna_mqc.tsv",
        sname=sample_name,
    )

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
