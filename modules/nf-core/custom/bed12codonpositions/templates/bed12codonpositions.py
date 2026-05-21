#!/usr/bin/env python3
"""Expand a BED12 file into a BED6 of in-frame positions along the spliced feature.

For each BED12 record, the script walks the block (exon) structure in
mRNA order and emits one BED6 row per in-frame position along the
spliced span, projecting each back to genomic coordinates. Frame is
defined relative to the start of the record itself.

BED12 convention (UCSC):
  chrom start end name score strand thickStart thickEnd itemRgb
  blockCount blockSizes blockStarts
where `start` is 0-based and `end` is half-open, blockStarts are offsets
from `start`, and blocks are listed in ascending genomic-coordinate
order on both strands. mRNA-order traversal is therefore left-to-right
on `+` strand and right-to-left on `-` strand.

Output (BED6, 0-based half-open):
  chrom  start  end  name  0  strand

By default the script emits the 5' nucleotide of every codon (step 3,
width 1). With `--width N` (N > 1) the script emits up to N consecutive
mRNA positions per codon; if those positions cross a block boundary the
positions are split into one BED row per block (so the union of rows
for a single codon still maps back to a contiguous mRNA region).
"""

import argparse
import sys


def expand_blocks(start, block_sizes, block_starts):
    blocks = [(start + off, start + off + sz) for sz, off in zip(block_sizes, block_starts)]
    blocks.sort()
    return blocks


def mrna_to_genomic_runs(blocks, strand, mrna_start, mrna_end):
    """Project a half-open mRNA span [mrna_start, mrna_end) onto genomic
    coordinates, returning a list of (g_start, g_end) BED-style runs
    (one per overlapped block, in ascending genomic order)."""
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

    runs.sort()
    return runs


def process(line, frame, step, width, keep_duplicates):
    parts = line.rstrip("\\n").split("\\t")
    if len(parts) < 12:
        return []
    chrom = parts[0]
    try:
        start = int(parts[1])
    except ValueError:
        return []
    name = parts[3]
    strand = parts[5]
    try:
        block_count = int(parts[9])
    except ValueError:
        return []
    block_sizes = [int(x) for x in parts[10].rstrip(",").split(",") if x != ""]
    block_starts = [int(x) for x in parts[11].rstrip(",").split(",") if x != ""]
    if len(block_sizes) != block_count or len(block_starts) != block_count:
        return []

    blocks = expand_blocks(start, block_sizes, block_starts)
    total_len = sum(be - bs for bs, be in blocks)

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
            rows.append(f"{chrom}\\t{g_start}\\t{g_end}\\t{name}\\t0\\t{strand}")
    return rows


def run(bed12_path, output_path, frame, step, width, keep_duplicates):
    with open(bed12_path) as fh, open(output_path, "w") as oh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            for row in process(line, frame, step, width, keep_duplicates):
                oh.write(row + "\\n")


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
    sys.stderr.write("--step must be positive\\n")
    sys.exit(2)
if parsed_args.width <= 0:
    sys.stderr.write("--width must be positive\\n")
    sys.exit(2)
if parsed_args.frame < 0:
    sys.stderr.write("--frame must be non-negative\\n")
    sys.exit(2)

run(
    "${bed12}",
    "${prefix}.bed",
    parsed_args.frame,
    parsed_args.step,
    parsed_args.width,
    parsed_args.keep_duplicates,
)

python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
with open("versions.yml", "w") as fh:
    fh.write('"${task.process}":\\n    python: "' + python_version + '"\\n')
