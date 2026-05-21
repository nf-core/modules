#!/usr/bin/env python3
"""Expand a BED12 into a BED6 of in-frame mRNA positions.

Walks each record's blocks in mRNA order (5'→3'), emits every --step-th
mRNA position starting at --frame, and projects them back to genomic
coordinates. Rows are written in mRNA-traversal order, so '-' strand
records come out in descending genomic order.
"""

import argparse
import platform
import sys

import pandas as pd
import yaml

BED12_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
]


def parse_block_field(value):
    return [int(x) for x in str(value).rstrip(",").split(",") if x != ""]


def mrna_to_genomic_runs(blocks, strand, mrna_start, mrna_end):
    """Project a half-open mRNA span [mrna_start, mrna_end) onto genomic
    coordinates, returning a list of (g_start, g_end) BED-style runs
    (one per overlapped block, in mRNA-traversal order)."""
    if strand == "+":
        ordered = list(blocks)
    elif strand == "-":
        ordered = list(reversed(blocks))
    else:
        return []

    runs = []
    cum = 0
    for blk_start, blk_end in ordered:
        blk_len = blk_end - blk_start
        blk_lo = cum
        blk_hi = cum + blk_len
        cum = blk_hi

        lo = max(mrna_start, blk_lo)
        hi = min(mrna_end, blk_hi)
        if lo >= hi:
            continue
        off_lo = lo - blk_lo
        off_hi = hi - blk_lo
        if strand == "+":
            g_lo = blk_start + off_lo
            g_hi = blk_start + off_hi
        else:
            g_hi = blk_end - off_lo
            g_lo = blk_end - off_hi
        runs.append((g_lo, g_hi))

    return runs


def emit_rows(row, frame, step, width, keep_duplicates):
    block_sizes = parse_block_field(row["blockSizes"])
    block_starts = parse_block_field(row["blockStarts"])
    if len(block_sizes) != int(row["blockCount"]) or len(block_starts) != int(row["blockCount"]):
        sys.stderr.write(
            f"warning: skipping {row['name']!r}: blockCount={row['blockCount']} but "
            f"blockSizes has {len(block_sizes)} entries and blockStarts has {len(block_starts)}\\n"
        )
        return []

    blocks = sorted((row["start"] + off, row["start"] + off + sz) for sz, off in zip(block_sizes, block_starts))
    total_len = sum(be - bs for bs, be in blocks)
    chrom = row["chrom"]
    name = row["name"]
    score = row["score"]
    strand = row["strand"]

    rows = []
    seen = set()
    for mrna_pos in range(frame, total_len, step):
        if mrna_pos + width > total_len:
            break
        for g_start, g_end in mrna_to_genomic_runs(blocks, strand, mrna_pos, mrna_pos + width):
            key = (chrom, g_start, g_end, name, strand)
            if not keep_duplicates and key in seen:
                continue
            seen.add(key)
            rows.append((chrom, g_start, g_end, name, score, strand))
    return rows


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "--frame",
    type=int,
    default=0,
    help="mRNA offset of the first position to emit (default: 0).",
)
parser.add_argument(
    "--step",
    type=int,
    default=3,
    help="Stride between successive emitted positions on the mRNA (default: 3).",
)
parser.add_argument(
    "--width",
    type=int,
    default=1,
    help="Width in nucleotides of each emitted span on the mRNA (default: 1). "
    "Spans that cross a block boundary are split into one BED row per block.",
)
parser.add_argument(
    "--keep-duplicates",
    action="store_true",
    help="Keep duplicate (chrom, start, end, name, strand) rows arising from "
    "the same record (e.g. when --width >= --step).",
)
parsed_args = parser.parse_args("${args}".split() if "${args}".strip() else [])

if parsed_args.step <= 0:
    raise SystemExit("--step must be positive")
if parsed_args.width <= 0:
    raise SystemExit("--width must be positive")
if parsed_args.frame < 0:
    raise SystemExit("--frame must be non-negative")

bed = pd.read_csv(
    "${bed12}",
    sep="\\t",
    comment="#",
    header=None,
    names=BED12_COLUMNS,
    dtype={"chrom": str, "name": str, "strand": str},
)
bed = bed[~bed["chrom"].astype(str).str.startswith(("track", "browser"))]

out_rows = []
for _, rec in bed.iterrows():
    out_rows.extend(
        emit_rows(
            rec,
            parsed_args.frame,
            parsed_args.step,
            parsed_args.width,
            parsed_args.keep_duplicates,
        )
    )

out = pd.DataFrame(out_rows, columns=["chrom", "start", "end", "name", "score", "strand"])
out.to_csv("${prefix}.bed", sep="\\t", header=False, index=False)

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
