#!/usr/bin/env python3
# Written by Jonathan Manning (@pinin4fjords). Released under the MIT license.

"""Collapse small ORFs sharing an amino-acid cluster into a single catalogue entry.

The coordinate-based merge in custom/orfmerge groups ORFs that overlap on the
genome, but the same micropeptide is frequently encoded at several distinct,
non-overlapping genomic loci (typically repetitive regions), and those copies
survive as separate catalogue rows. Following the GENCODE Ribo-seq ORF catalogue
convention (Mudge et al. 2022, Nat Biotechnol, doi:10.1038/s41587-022-01369-0;
gencode-riboseqORFs collapse_cutoff 0.9), small ORFs (orf_class == "smORF", i.e.
aa_length <= 100) are clustered by amino-acid sequence identity upstream
(mmseqs/easycluster) and this module folds each multi-member cluster down to one
representative.

Only smORF rows are collapsed; larger ORFs and transcript-anchored classes pass
through untouched, preserving the deterministic coordinate/transcript merge from
upstream. Among the smORF members of a cluster the representative is chosen here
(longest aa_length, ties broken by orf_id) so the result is independent of which
sequence MMseqs2 labelled the cluster representative. Catalogue row order is
preserved; dropped members fold their cross-caller / cross-sample evidence and
gene mappings into the survivor.
"""

import argparse
import platform
import re
import shlex
import sys
from collections import defaultdict

import pandas as pd
import yaml

CALLERS = ("ribotish", "ribocode", "ribotricer", "rpbp", "price")
# `bedtools getfasta -nameOnly -s` appends the strand as "(+)"/"(-)" to each
# sequence name, so FASTA headers and MMseqs2 cluster ids arrive as
# "<orf_id>(+)". Strip it to recover the bare orf_id used in the catalogue.
STRAND_SUFFIX = re.compile(r"\\([+-]\\)\$")
SCORE_DIRECTIONS = {
    "ribotish": "min",
    "ribocode": "min",
    "ribotricer": "max",
    "rpbp": "max",
    "price": "min",
}
CLASS_ORDER = ("canonical_cds", "uORF", "dORF", "novel_u", "smORF", "other")
SMORF_CLASS = "smORF"


def read_fasta(path):
    seqs = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                header = line[1:].split()[0] if len(line) > 1 else ""
                name = STRAND_SUFFIX.sub("", header)
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_clusters(path):
    """Map member orf_id -> cluster key (representative id) from an mmseqs cluster TSV."""
    cluster_of = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) >= 2:
                member = STRAND_SUFFIX.sub("", fields[1])
                rep = STRAND_SUFFIX.sub("", fields[0])
                cluster_of[member] = rep
    return cluster_of


def best_score(values, direction):
    nums = []
    for v in values:
        try:
            if v not in ("", None):
                nums.append(float(v))
        except ValueError:
            continue
    if not nums:
        return ""
    return f"{(max(nums) if direction == 'max' else min(nums)):.6g}"


def merge_members(members):
    """Fold smORF rows sharing an AA cluster into one representative row dict."""
    rep = sorted(members, key=lambda r: (-int(r.get("aa_length") or 0), r["orf_id"]))[0]
    out = dict(rep)
    for c in CALLERS:
        out[f"called_by_{c}"] = (
            "1" if any(r.get(f"called_by_{c}") == "1" for r in members) else "0"
        )
        out[f"score_{c}"] = best_score(
            [r.get(f"score_{c}", "") for r in members], SCORE_DIRECTIONS[c]
        )
    samples = sorted(
        {s for r in members for s in (r.get("samples") or "").split(",") if s}
    )
    out["n_samples"] = str(len(samples))
    out["samples"] = ",".join(samples)
    return rep["orf_id"], out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(shlex.split("${args}"))

    prefix = "${prefix}"

    catalogue = pd.read_csv(
        "${catalogue_tsv}", sep="\\t", comment="#", dtype=str, keep_default_na=False
    )
    header = list(catalogue.columns)
    rows = catalogue.to_dict("records")

    bed_index = {}
    with open("${bed12}") as fh:
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) >= 12:
                bed_index[parts[3]] = line.rstrip("\\n")

    aa = read_fasta("${aa_fasta}")
    cluster_of = read_clusters("${cluster_tsv}")

    clusters = defaultdict(list)
    for r in rows:
        if r.get("orf_class") == SMORF_CLASS:
            clusters[cluster_of.get(r["orf_id"], r["orf_id"])].append(r)

    remap, merged_rows, dropped = {}, {}, set()
    for members in clusters.values():
        if len(members) < 2:
            continue
        rep_id, merged = merge_members(members)
        merged_rows[rep_id] = merged
        for m in members:
            remap[m["orf_id"]] = rep_id
            if m["orf_id"] != rep_id:
                dropped.add(m["orf_id"])

    kept = [merged_rows.get(r["orf_id"], r) for r in rows if r["orf_id"] not in dropped]
    out = pd.DataFrame(kept, columns=header)
    out.to_csv(f"{prefix}.tsv", sep="\\t", index=False)

    with (
        open(f"{prefix}.bed12", "w") as bh,
        open(f"{prefix}.fasta", "w") as ah,
    ):
        for oid in out["orf_id"]:
            if oid in bed_index:
                bh.write(bed_index[oid] + "\\n")
            if oid in aa:
                ah.write(f">{oid}\\n{aa[oid]}\\n")

    o2g = pd.read_csv("${orf_to_gene_tsv}", sep="\\t", dtype=str, keep_default_na=False)
    o2g["orf_id"] = o2g["orf_id"].map(lambda o: remap.get(o, o))
    o2g.drop_duplicates().to_csv(f"{prefix}.orf_to_gene.tsv", sep="\\t", index=False)

    counts = out["orf_class"].value_counts()
    with open(f"{prefix}.mqc.tsv", "w") as mh:
        mh.write(
            "# id: orf_catalogue\\n"
            "# section_name: 'ORF catalogue'\\n"
            "# description: 'Per-class ORF counts in the merged catalogue.'\\n"
            "# plot_type: 'table'\\n"
            "# pconfig:\\n"
            "#   id: 'orf_catalogue_table'\\n"
            "#   title: 'ORF catalogue'\\n"
            "Class\\tCount\\n"
        )
        for cls in CLASS_ORDER:
            mh.write(f"{cls}\\t{int(counts.get(cls, 0))}\\n")

    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {
                "${task.process}": {
                    "python": platform.python_version(),
                    "pandas": pd.__version__,
                }
            },
            fh,
            default_flow_style=False,
            sort_keys=False,
        )
    return 0


sys.exit(main())
