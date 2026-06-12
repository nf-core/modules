#!/usr/bin/env python3
"""Merge per-sample, per-caller normalised ORF BED12s into a unified catalogue.

Strategy is class-aware (operating on the harmonised `orf_class` from
custom/orfnormalise):

  - canonical_cds, uORF, dORF, other:  collapse by transcript_id (handles
    multi-exon CDS correctly; uORF/dORF are anchored to a host transcript).
  - novel_u, smORF: greedy reciprocal-overlap clustering on the outer
    genomic span at `--reciprocal-overlap` (default 0.8). Catches fuzzy
    cross-caller matches and exact-coordinate collapses in one pass.

Cross-caller consensus is recorded in two column families on the output
catalogue TSV:

  - `called_by_<caller>`: 0/1 indicator per supported caller.
  - `score_<caller>`:     best score from that caller within the cluster.
                          Score direction is per-caller (p-values are
                          minimised; Bayes factors / phase scores are
                          maximised).

Cross-sample recurrence is recorded in two further columns:

  - `n_samples`: number of distinct samples contributing to the cluster.
  - `samples`:   sorted, comma-separated list of those sample ids.

Outputs:

  ${prefix}.catalogue.bed12      merged catalogue (genomic blocks).
  ${prefix}.catalogue.tsv        per-ORF table with caller-tracking cols.
  ${prefix}.orf_to_gene.tsv      one row per (orf_id, gene_id, transcript_id);
                                 an ORF can map to multiple host transcripts.
  ${prefix}.catalogue.mqc.tsv    MultiQC custom-content sidecar
                                 (per-class counts).
"""

import argparse
import csv
import glob
import platform
import shlex
import sys
from collections import defaultdict
from pathlib import Path

import yaml

CALLERS = ("ribotish", "ribocode", "ribotricer", "rpbp", "price")

# Direction of each caller's native score: 'min' = lower is better
# (p-values); 'max' = higher is better (Bayes factors, phase scores).
# When aggregating per-cluster scores in the catalogue TSV, we keep the
# direction-appropriate best.
SCORE_DIRECTIONS = {
    "ribotish": "min",
    "ribocode": "min",
    "ribotricer": "max",
    "rpbp": "max",
    "price": "min",
}

CLASS_ORDER = ("canonical_cds", "uORF", "dORF", "novel_u", "smORF", "other")


def group_by(rows, keyfn):
    """Group rows into clusters keyed by ``keyfn(row)``, preserving first-seen order."""
    groups = defaultdict(list)
    for r in rows:
        groups[keyfn(r)].append(r)
    return list(groups.values())


def cluster_by_reciprocal_overlap(rows, frac=0.8):
    """Greedy clustering by reciprocal overlap >= ``frac`` on the outer span.

    Used to merge cross-caller calls of the same biological ORF when they
    don't share a transcript_id (e.g. novel intergenic from different
    callers). O(N^2) worst case but bounded by per-run ORF counts.

    Greedy and order-dependent: at the default frac=0.8 a chain A-B-C
    where A overlaps B at 0.85, B overlaps C at 0.85, but A overlaps C
    at 0.75 can either land as {A, B, C} or {A, B} + {C} depending on
    iteration order. Acceptable in practice (overlap chains that
    straddle the threshold are rare at frac=0.8) but worth flagging for
    consumers that want strictly transitive clustering.
    """
    clusters = []
    assigned = [False] * len(rows)
    order = sorted(
        range(len(rows)),
        key=lambda i: (rows[i]["chrom"], rows[i]["strand"], int(rows[i]["start"])),
    )
    for i in order:
        if assigned[i]:
            continue
        ri = rows[i]
        ri_start, ri_end = int(ri["start"]), int(ri["end"])
        ri_len = ri_end - ri_start
        cluster = [ri]
        assigned[i] = True
        for j in order:
            if assigned[j]:
                continue
            rj = rows[j]
            if rj["chrom"] != ri["chrom"] or rj["strand"] != ri["strand"]:
                continue
            rj_start, rj_end = int(rj["start"]), int(rj["end"])
            if rj_start >= ri_end or rj_end <= ri_start:
                continue
            ov = min(ri_end, rj_end) - max(ri_start, rj_start)
            if ov <= 0:
                continue
            rj_len = rj_end - rj_start
            if ri_len > 0 and rj_len > 0 and ov / ri_len >= frac and ov / rj_len >= frac:
                cluster.append(rj)
                assigned[j] = True
        clusters.append(cluster)
    return clusters


def representative(cluster):
    """Pick a representative row from a cluster.

    Preference order: canonical_cds, then uORF/dORF, then novel_u/smORF,
    then other; ties broken by longest aa_length.
    """
    rank = {"canonical_cds": 0, "uORF": 1, "dORF": 1, "novel_u": 2, "smORF": 2, "other": 3}
    return sorted(
        cluster,
        key=lambda r: (rank.get(r.get("orf_class", "other"), 3), -int(r.get("aa_length") or 0)),
    )[0]


def aggregate_caller_columns(cluster):
    out = {}
    by_caller = defaultdict(list)
    for r in cluster:
        by_caller[r.get("caller", "")].append(r)
    for c in CALLERS:
        out[f"called_by_{c}"] = "1" if by_caller.get(c) else "0"
    for c in CALLERS:
        rows_c = by_caller.get(c, [])
        if not rows_c:
            out[f"score_{c}"] = ""
            continue
        scores = []
        for r in rows_c:
            s = r.get("score", "")
            try:
                if s != "":
                    scores.append(float(s))
            except ValueError:
                continue
        if not scores:
            out[f"score_{c}"] = ""
        else:
            best = max(scores) if SCORE_DIRECTIONS[c] == "max" else min(scores)
            out[f"score_{c}"] = f"{best:.6g}"
    return out


def load_normalised(tsv_paths, bed_paths):
    rows = []
    bed_index = {}
    for p in tsv_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p, newline="") as fh:
            data = (line for line in fh if not line.startswith("#"))
            rows.extend(csv.DictReader(data, delimiter="\\t"))
    for p in bed_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        with open(p) as fh:
            for line in fh:
                parts = line.rstrip("\\n").split("\\t")
                if len(parts) < 12:
                    continue
                orf_id = parts[3]
                if orf_id not in bed_index:
                    bed_index[orf_id] = line.rstrip("\\n")
    return rows, bed_index


def write_catalogue(prefix, clusters, bed_index):
    cat_bed = Path(f"{prefix}.catalogue.bed12")
    cat_tsv = Path(f"{prefix}.catalogue.tsv")
    o2g_tsv = Path(f"{prefix}.orf_to_gene.tsv")
    mqc_tsv = Path(f"{prefix}.catalogue.mqc.tsv")

    catalogue_cols = (
        ["orf_id", "chrom", "start", "end", "strand", "gene_id", "transcript_id", "orf_class", "aa_length"]
        + [f"called_by_{c}" for c in CALLERS]
        + [f"score_{c}" for c in CALLERS]
        + ["n_samples", "samples"]
    )

    per_class_counts = defaultdict(int)

    with open(cat_bed, "w") as bh, open(cat_tsv, "w") as th, open(o2g_tsv, "w") as oh:
        th.write("\\t".join(catalogue_cols) + "\\n")
        oh.write("orf_id\\tgene_id\\ttranscript_id\\n")
        for idx, cluster in enumerate(clusters):
            rep = representative(cluster)
            stable_id = f"orf_{idx + 1:08d}"
            bed_line = bed_index.get(rep["orf_id"])
            if bed_line:
                parts = bed_line.split("\\t")
                parts[3] = stable_id
                bh.write("\\t".join(parts) + "\\n")

            caller_cols = aggregate_caller_columns(cluster)
            row_out = [
                stable_id,
                rep.get("chrom", ""),
                rep.get("start", ""),
                rep.get("end", ""),
                rep.get("strand", ""),
                rep.get("gene_id", ""),
                rep.get("transcript_id", ""),
                rep.get("orf_class", "other"),
                rep.get("aa_length", "0"),
            ]
            row_out += [caller_cols[f"called_by_{c}"] for c in CALLERS]
            row_out += [caller_cols[f"score_{c}"] for c in CALLERS]

            sample_ids = sorted({r.get("sample_id", "") for r in cluster if r.get("sample_id")})
            row_out += [str(len(sample_ids)), ",".join(sample_ids)]
            th.write("\\t".join(row_out) + "\\n")

            per_class_counts[rep.get("orf_class", "other")] += 1

            seen_gt = set()
            for r in cluster:
                key = (r.get("gene_id", ""), r.get("transcript_id", ""))
                if key in seen_gt:
                    continue
                seen_gt.add(key)
                oh.write(f"{stable_id}\\t{key[0]}\\t{key[1]}\\n")

    with open(mqc_tsv, "w") as mh:
        mh.write("# id: orf_catalogue\\n")
        mh.write("# section_name: 'ORF catalogue'\\n")
        mh.write("# description: 'Per-class ORF counts in the merged catalogue.'\\n")
        mh.write("# plot_type: 'table'\\n")
        mh.write("# pconfig:\\n")
        mh.write("#   id: 'orf_catalogue_table'\\n")
        mh.write("#   title: 'ORF catalogue'\\n")
        mh.write("Class\\tCount\\n")
        for cls in CLASS_ORDER:
            mh.write(f"{cls}\\t{per_class_counts.get(cls, 0)}\\n")


def write_versions():
    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {"${task.process}": {"python": platform.python_version()}},
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        default=0.8,
        help="Reciprocal-overlap fraction for novel_u/smORF clustering (default: 0.8)",
    )
    args = parser.parse_args(shlex.split("${args}"))

    bed_paths = sorted(Path(p) for p in glob.glob("beds/*"))
    tsv_paths = sorted(Path(p) for p in glob.glob("tsvs/*"))

    prefix = "${prefix}"
    rows, bed_index = load_normalised(tsv_paths, bed_paths)

    if not rows:
        write_catalogue(prefix, [], {})
        write_versions()
        return 0

    by_class = defaultdict(list)
    for r in rows:
        by_class[r.get("orf_class", "other")].append(r)

    clusters = []
    # canonical CDS: one per transcript by definition - collapse by (tid, strand).
    clusters.extend(group_by(by_class.get("canonical_cds", []), lambda r: (r.get("transcript_id") or "", r["strand"])))
    # uORF/dORF/other: a transcript can host multiple distinct ones, so
    # additionally key on the outer span to keep them separate.
    for cls in ("uORF", "dORF", "other"):
        clusters.extend(
            group_by(
                by_class.get(cls, []),
                lambda r: (r.get("transcript_id") or "", r["strand"], int(r["start"]), int(r["end"])),
            )
        )
    # novel_u / smORF: not transcript-anchored - reciprocal-overlap clustering.
    for cls in ("novel_u", "smORF"):
        clusters.extend(cluster_by_reciprocal_overlap(by_class.get(cls, []), frac=args.reciprocal_overlap))

    write_catalogue(prefix, clusters, bed_index)
    write_versions()
    return 0


sys.exit(main())
