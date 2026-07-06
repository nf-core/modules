#!/usr/bin/env python3
# Written by Jonathan Manning (@pinin4fjords). Released under the MIT license.

"""Collapse small ORFs sharing an amino-acid cluster into a single catalogue entry.

The coordinate-based merge in custom/orfmerge groups ORFs that overlap on the
genome, but the same micropeptide is frequently encoded at several distinct,
non-overlapping genomic loci (typically repetitive regions), and those copies
survive as separate catalogue rows. This peptide-level deduplication is inspired
by the GENCODE Ribo-seq ORF consolidation (Mudge et al. 2022, Nat Biotechnol,
doi:10.1038/s41587-022-01369-0; gencode-riboseqORFs collapse_cutoff 0.9), with
two deliberate departures from that reference:

  - GENCODE collapses overlapping ORFs of any size within a shared locus;
    this restricts collapsing to small ORFs (orf_class == "smORF", i.e.
    aa_length <= 100) and clusters them locus-agnostically across the whole
    catalogue, since the target case is one micropeptide recurring at several
    non-overlapping loci. The smORF-only restriction is this pipeline's choice,
    not a GENCODE property.
  - similarity is MMseqs2 global sequence identity (--min-seq-id 0.9,
    mmseqs/easycluster upstream) rather than GENCODE's longest-shared-substring
    / P-site-overlap metric, so the 0.9 here approximates rather than
    reproduces collapse_cutoff 0.9.

Small ORFs are clustered by amino-acid identity upstream (mmseqs/easycluster)
and this module folds each multi-member cluster down to one representative.

Alongside the collapsed catalogue it emits a consensus view (`*.consensus.*`)
filtered to ORFs supported by at least `--min-callers` distinct callers and
recurring in at least `--min-samples` samples (both default 1, i.e. no
filtering). The filter is applied after the collapse, so the consensus is the
high-confidence subset of the de-redundified catalogue and a folded
micropeptide is judged on its combined cross-caller / cross-sample evidence.

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
        out[f"called_by_{c}"] = "1" if any(r.get(f"called_by_{c}") == "1" for r in members) else "0"
        out[f"score_{c}"] = best_score([r.get(f"score_{c}", "") for r in members], SCORE_DIRECTIONS[c])
    samples = sorted({s for r in members for s in (r.get("samples") or "").split(",") if s})
    out["n_samples"] = str(len(samples))
    out["samples"] = ",".join(samples)
    return rep["orf_id"], out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--min-callers",
        type=int,
        default=1,
        help="Minimum distinct callers for an ORF to enter the consensus view (default: 1)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum distinct samples for an ORF to enter the consensus view (default: 1)",
    )
    args = parser.parse_args(shlex.split("${args}"))

    prefix = "${prefix}"

    catalogue = pd.read_csv("${catalogue_tsv}", sep="\\t", comment="#", dtype=str, keep_default_na=False)
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

    o2g = pd.read_csv("${orf_to_gene_tsv}", sep="\\t", dtype=str, keep_default_na=False)
    o2g["orf_id"] = o2g["orf_id"].map(lambda o: remap.get(o, o))
    o2g = o2g.drop_duplicates()

    def write_view(out_prefix, df):
        """Write the catalogue TSV, BED12 and ORF-to-gene mapping for `df`."""
        df.to_csv(f"{out_prefix}.tsv", sep="\\t", index=False)
        with open(f"{out_prefix}.bed12", "w") as bh:
            for oid in df["orf_id"]:
                if oid in bed_index:
                    bh.write(bed_index[oid] + "\\n")
        o2g[o2g["orf_id"].isin(set(df["orf_id"]))].to_csv(f"{out_prefix}.orf_to_gene.tsv", sep="\\t", index=False)

    # Full collapsed catalogue (+ amino-acid FASTA).
    write_view(prefix, out)
    with open(f"{prefix}.fasta", "w") as ah:
        for oid in out["orf_id"]:
            if oid in aa:
                ah.write(f">{oid}\\n{aa[oid]}\\n")

    # Consensus view: high-confidence subset of the collapsed catalogue, so a
    # folded micropeptide is judged on its combined cross-caller / cross-sample
    # evidence. At the default thresholds (1, 1) it equals the full catalogue.
    n_callers = sum((out[f"called_by_{c}"] == "1").astype(int) for c in CALLERS)
    n_samples = pd.to_numeric(out["n_samples"], errors="coerce").fillna(0).astype(int)
    consensus = out[(n_callers >= args.min_callers) & (n_samples >= args.min_samples)]
    write_view(f"{prefix}.consensus", consensus)

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
                    "pyyaml": yaml.__version__,
                }
            },
            fh,
            default_flow_style=False,
            sort_keys=False,
        )
    return 0


sys.exit(main())
