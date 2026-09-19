#!/usr/bin/env python3
"""Normalise an ORF caller's per-sample output to a unified BED12 + sidecar TSV.

One template, five callers. The `caller` val input (one of ribocode,
ribotish, ribotricer, rpbp, price) selects the per-caller parser,
classifier and score mapper.

`orf_class` is purely positional: it records where the ORF sits relative to
the annotated CDS and never encodes its length. Use the `is_smorf` column to
select small ORFs.

Harmonised orf_class vocabulary written into the sidecar TSV:
    canonical_cds: ORF maps to an annotated CDS (including truncated /
                   extended variants of one).
    uORF:          upstream ORF, not overlapping the CDS.
    uoORF:         upstream ORF overlapping the CDS out of frame.
    dORF:          downstream ORF, not overlapping the CDS.
    doORF:         downstream ORF overlapping the CDS out of frame.
    intORF:        ORF contained within the CDS, out of frame.
    novel_u:       novel / intergenic ORF not assigned to an annotated CDS.
    other:         anything a caller cannot place, including ribotricer's
                   `internal` fall-through and PRICE's `intronic`.

Not every caller can report every class. Ribo-TISH's `5'UTR` covers both uORF
and uoORF (tisType() tests only the start position) and PRICE's enum has no
doORF, so absence of a class for one caller is a tool limitation rather than
evidence about the ORF. `orf_type_native` carries the caller's own label so
these decisions stay auditable.

Source-column choice for the per-ORF score, ORF-type classification, and
length-derivation is configurable via ext.args:

    --score-field NAME       Override the column read into `score`.
    --orf-type-field NAME    Override the column read into `orf_type`
                             (before classification).
    --length-field NAME      Override the nucleotide-length column used to
                             derive `aa_length` (applies to callers that
                             derive length from a column: ribocode,
                             ribotricer, rpbp).
    --aa-length-field NAME   Override a direct aa-length column (ribotish
                             only - the only caller emitting AALen).
    --smorf-max-aa N         Maximum aa_length flagged in `is_smorf`
                             (default 100). Affects the flag only, never
                             `orf_class`.

Defaults are per-caller and listed in DEFAULT_FIELDS below; if a column
chain is provided, the parser walks it in order. The resolved column
name for each field is written into the TSV header as a
`# parser_columns: ...` comment line for downstream provenance.
"""

import argparse
import csv
import gzip
import io
import math
import platform
import re
import shlex
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import yaml

csv.field_size_limit(sys.maxsize)

# Bumped whenever the orf_class vocabulary changes, so a published catalogue
# self-identifies which vocabulary produced it.
ORF_CLASS_VOCABULARY = 2

CALLER = "${caller}"
INPUT = Path("${orfs_table}")
GTF = Path("${gtf}")
SAMPLE_ID = "${sample_id}"
OUT_BED = Path("${prefix}.bed12")
OUT_TSV = Path("${prefix}.tsv")

TSV_HEADER = (
    "orf_id\\tcaller\\tsample_id\\tchrom\\tstart\\tend\\tstrand\\t"
    "gene_id\\ttranscript_id\\torf_class\\taa_length\\tscore\\torf_type_native\\tis_smorf"
)

# Per-caller default field-name preference chains. Each list is walked in
# order until a row supplies a usable (non-empty, non-"None") value.
# `None` means the field is not column-read for that caller (derived
# elsewhere or absent from the source format).
DEFAULT_FIELDS = {
    "ribocode": {
        "score": ["pval_combined", "Pval_combined"],
        "orf_type": ["ORF_type", "Type"],
        "length": ["ORF_length", "ORFlength"],
        "aa_length": None,
    },
    "ribotish": {
        "score": ["FisherPvalue", "Pvalcombined", "RiboPvalue", "TISPvalue", "Pvalue"],
        "orf_type": ["TisType", "TISType"],
        "length": None,
        "aa_length": ["AALen"],
    },
    "ribotricer": {
        "score": ["phase_score"],
        "orf_type": ["ORF_type"],
        "length": ["length"],
        "aa_length": None,
    },
    "rpbp": {
        "score": ["bayes_factor_mean"],
        "orf_type": None,
        "length": ["orf_len"],
        "aa_length": None,
    },
    "price": {
        "score": ["p value", "p_value"],
        "orf_type": ["Type"],
        "length": None,
        "aa_length": None,
    },
}

# RPBP predicted-orfs BED column names (the file ships a `#`-prefixed
# header but we keep our internal names clean).
RPBP_COLUMNS = [
    "seqname",
    "start",
    "end",
    "id",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "color",
    "num_exons",
    "exon_lengths",
    "exon_genomic_relative_starts",
    "orf_num",
    "orf_len",
    "p_translated_mean",
    "p_translated_var",
    "p_background_mean",
    "p_background_var",
    "bayes_factor_mean",
    "bayes_factor_var",
    "chi_square_p",
    "x_1_sum",
    "x_2_sum",
    "x_3_sum",
    "profile_sum",
]


_ATTR_RE = re.compile(r'(\\w+)\\s+"([^"]*)"')


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: list

    @property
    def length(self):
        return sum(e[1] - e[0] for e in self.exons)


# ----------------------------------------------------------------------------
# Small shared helpers
# ----------------------------------------------------------------------------


def clamp1000(x):
    return max(0, min(1000, int(round(x))))


def _resolve_chain(caller, field_key, override):
    if override:
        return [override]
    return DEFAULT_FIELDS[caller].get(field_key) or []


def pick(row, fields, resolved, key):
    """Return the first usable value for `key` from this row, walking the
    resolved column chain in `fields[key]`. Records the column actually read
    in `resolved[key]` for provenance. Returns None if nothing usable."""
    for col in fields.get(key) or []:
        val = row.get(col)
        if val is None or val == "" or val == "None":
            continue
        resolved[key] = col
        return val
    return None


def aa_from_length(row, fields, resolved):
    """Derive aa_length from a nucleotide-length column (ribocode/ribotricer/rpbp)."""
    raw = pick(row, fields, resolved, "length")
    try:
        nt = int(raw) if raw else 0
    except (TypeError, ValueError):
        nt = 0
    return max(0, (nt - 3) // 3) if nt > 0 else 0


def parse_intervals(s, seps=("-",)):
    """Parse comma-separated 1-based inclusive `start<sep>end` tokens into
    sorted 0-based half-open (start, end) blocks. `seps` is tried in order
    per token (ribotish blocks use ':' or '-'; ribotricer coords use '-')."""
    if not s:
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        pair = None
        for sep in seps:
            if sep in tok:
                pair = tok.split(sep, 1)
                break
        if pair is None:
            continue
        try:
            a_i, b_i = int(pair[0]), int(pair[1])
        except ValueError:
            continue
        if b_i < a_i:
            a_i, b_i = b_i, a_i
        out.append((a_i - 1, b_i))
    out.sort()
    return out


# ----------------------------------------------------------------------------
# I/O + GTF helpers
# ----------------------------------------------------------------------------


def open_text(path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    return open(p, encoding="utf-8")


def parse_attributes(field):
    return dict(_ATTR_RE.findall(field))


def load_transcripts(gtf_path):
    p = Path(gtf_path)
    if not p.is_file() or p.stat().st_size == 0:
        return {}
    by_tid = {}
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom = parts[0]
            start_1 = int(parts[3])
            end_1 = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            tx = by_tid.get(tid)
            if tx is None:
                tx = Transcript(
                    transcript_id=tid,
                    gene_id=attrs.get("gene_id", ""),
                    chrom=chrom,
                    strand=strand,
                    exons=[],
                )
                by_tid[tid] = tx
            tx.exons.append((start_1 - 1, end_1))
    for tx in by_tid.values():
        tx.exons.sort()
    return by_tid


def transcript_to_genomic_blocks(tx, tx_start_0, tx_end_excl):
    if tx_end_excl <= tx_start_0:
        return []
    length = tx.length
    if tx_start_0 < 0 or tx_end_excl > length:
        tx_start_0 = max(0, tx_start_0)
        tx_end_excl = min(length, tx_end_excl)
        if tx_end_excl <= tx_start_0:
            return []
    blocks = []
    if tx.strand == "+":
        cumulative = 0
        for gs, ge in tx.exons:
            esize = ge - gs
            exon_mrna_start = cumulative
            exon_mrna_end = cumulative + esize
            lo = max(tx_start_0, exon_mrna_start)
            hi = min(tx_end_excl, exon_mrna_end)
            if hi > lo:
                gstart = gs + (lo - exon_mrna_start)
                gend = gs + (hi - exon_mrna_start)
                blocks.append((gstart, gend))
            cumulative = exon_mrna_end
    else:
        cumulative = 0
        for gs, ge in reversed(tx.exons):
            esize = ge - gs
            exon_mrna_start = cumulative
            exon_mrna_end = cumulative + esize
            lo = max(tx_start_0, exon_mrna_start)
            hi = min(tx_end_excl, exon_mrna_end)
            if hi > lo:
                gend = ge - (lo - exon_mrna_start)
                gstart = ge - (hi - exon_mrna_start)
                blocks.append((gstart, gend))
            cumulative = exon_mrna_end
        blocks.sort()
    return blocks


def emit_bed12(chrom, blocks, name, score, strand):
    if not blocks:
        return ""
    start = blocks[0][0]
    end = blocks[-1][1]
    block_count = len(blocks)
    block_sizes = ",".join(str(b[1] - b[0]) for b in blocks) + ","
    block_starts = ",".join(str(b[0] - start) for b in blocks) + ","
    return "\\t".join(
        [
            chrom,
            str(start),
            str(end),
            name,
            str(score),
            strand,
            str(start),
            str(end),
            "0",
            str(block_count),
            block_sizes,
            block_starts,
        ]
    )


def emit_tsv_row(
    orf_id,
    caller,
    sample_id,
    chrom,
    start,
    end,
    strand,
    gene_id,
    transcript_id,
    orf_class,
    aa_length,
    score,
    orf_type_native,
    is_smorf,
):
    return "\\t".join(
        [
            orf_id,
            caller,
            sample_id,
            chrom,
            str(start),
            str(end),
            strand,
            gene_id,
            transcript_id,
            orf_class,
            str(aa_length),
            str(score),
            orf_type_native,
            "1" if is_smorf else "0",
        ]
    )


def write_outputs(bed_path, tsv_path, bed_lines, tsv_rows, parser_columns, unmapped, smorf_max_aa):
    bed_lines = [line for line in bed_lines if line]
    with open(bed_path, "w") as bh:
        for line in bed_lines:
            bh.write(line + "\\n")
    with open(tsv_path, "w") as th:
        pc = " ".join(f"{k}={v}" for k, v in parser_columns.items() if v)
        # An unmapped count above zero means this caller emitted an ORF-type
        # token CLASS_TOKENS does not know, so those rows fell into `other`.
        provenance = (
            f"# parser_columns: caller={CALLER} {pc}"
            f" unmapped_orf_type={sum(unmapped.values())}"
            f" orf_class_vocabulary={ORF_CLASS_VOCABULARY} smorf_max_aa={smorf_max_aa}"
        )
        if unmapped:
            provenance += f" unmapped_tokens={','.join(sorted(unmapped))}"
        th.write(provenance + "\\n")
        th.write(TSV_HEADER + "\\n")
        for r in tsv_rows:
            th.write(r + "\\n")


# ----------------------------------------------------------------------------
# ORF-type classification
#
# Every caller's ORF-type vocabulary is a closed enum, so tokens are matched
# exactly (casefolded) rather than by substring. Substring matching mis-fired
# on the overlap forms, because "uorf" is a substring of "overlap_uorf".
#
# No key may contain a colon: classify() matches only the part before one.
# ----------------------------------------------------------------------------

CLASS_TOKENS = {
    # RiboCode detectORF.classfy_orf
    "ribocode": {
        "annotated": "canonical_cds",
        "ccds": "canonical_cds",
        "uorf": "uORF",
        "dorf": "dORF",
        "overlap_uorf": "uoORF",
        "overlap_dorf": "doORF",
        "internal": "intORF",
        "novel": "novel_u",
        "intergenic": "novel_u",
    },
    # Ribo-TISH predict.TIS_types. 5'UTR conflates uORF with uoORF and
    # 3'UTR conflates dORF with doORF; tisType() only tests the start
    # position, so the overlap forms are not recoverable here.
    "ribotish": {
        "annotated": "canonical_cds",
        "truncated": "canonical_cds",
        "extended": "canonical_cds",
        "5'utr": "uORF",
        "uorf": "uORF",
        "3'utr": "dORF",
        "dorf": "dORF",
        "internal": "intORF",
        "novel": "novel_u",
        "intergenic": "novel_u",
    },
    # ribotricer prepare_orfs.check_orf_type. `internal` is that function's
    # terminal fall-through rather than a positive out-of-frame call, so it
    # stays `other`: mapping it to intORF would assert a frame relationship
    # ribotricer never tested.
    "ribotricer": {
        "annotated": "canonical_cds",
        "ccds": "canonical_cds",
        "uorf": "uORF",
        "super_uorf": "uORF",
        "overlap_uorf": "uoORF",
        "dorf": "dORF",
        "super_dorf": "dORF",
        "overlap_dorf": "doORF",
        "internal": "other",
        "novel": "novel_u",
        "intergenic": "novel_u",
    },
    # Rp-Bp defaults.orf_type_name_map. Only "canonical" is reachable today
    # because parse_rpbp has no type column to read; the rest of the
    # vocabulary arrives with the orfs-labels file.
    "rpbp": {
        "canonical": "canonical_cds",
        "annotated": "canonical_cds",
        "five_prime": "uORF",
        "uorf": "uORF",
        "five_prime_overlap": "uoORF",
        "three_prime": "dORF",
        "dorf": "dORF",
        "three_prime_overlap": "doORF",
        "internal": "intORF",
        "novel": "novel_u",
        "intergenic": "novel_u",
    },
    # PRICE PriceOrfType. The enum has no doORF. `orphan` is both the
    # not-transcript-consistent label and the final fall-through, so it maps
    # to novel_u; `intronic` has no positional equivalent and stays `other`.
    "price": {
        "cds": "canonical_cds",
        "ext": "canonical_cds",
        "trunc": "canonical_cds",
        "variant": "canonical_cds",
        "uorf": "uORF",
        "uoorf": "uoORF",
        "dorf": "dORF",
        "iorf": "intORF",
        "ncrna": "novel_u",
        "intronic": "other",
        "orphan": "novel_u",
    },
}


def classify(caller, orf_type):
    """Map a caller's native ORF-type token to the harmonised class.

    Only the part before the first colon is matched; the full token survives in
    `orf_type_native`.

    Returns (orf_class, matched); `matched` is False when a non-empty token's
    location matched no entry, so unmapped labels can be counted rather than
    silently absorbed into `other`.
    """
    if not orf_type:
        return "other", True
    location = orf_type.split(":", 1)[0]
    cls = CLASS_TOKENS[caller].get(location.strip().casefold())
    if cls is None:
        return "other", False
    return cls, True


# ----------------------------------------------------------------------------
# Per-caller parsers
#
# Each parser receives (path, transcripts, fields) where `fields` is a dict
# of resolved column-name chains for this caller. Returns a list of ORF dicts
# plus a `resolved` dict recording the column actually read for each field
# (for provenance reporting).
# ----------------------------------------------------------------------------


def parse_ribocode(path, transcripts, fields):
    resolved = {"score": "", "orf_type": "", "length": ""}
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows, resolved
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return rows, resolved
        for row in reader:
            tid = row.get("transcript_id") or row.get("Transcript_id") or ""
            gid = row.get("gene_id") or row.get("Gene_id") or ""
            chrom = row.get("chrom") or row.get("Chrom") or ""
            strand = row.get("strand") or row.get("Strand") or "+"
            try:
                t_start_1 = int(row.get("ORF_tstart") or row.get("ORF_Tstart") or 0)
                t_stop_1 = int(row.get("ORF_tstop") or row.get("ORF_Tstop") or 0)
            except ValueError:
                continue

            aa_len = aa_from_length(row, fields, resolved)
            orf_type = pick(row, fields, resolved, "orf_type") or ""

            score_raw = pick(row, fields, resolved, "score")
            try:
                pval = float(score_raw) if score_raw else 1.0
                bed_score = clamp1000((1.0 - pval) * 1000)
            except (TypeError, ValueError):
                pval = float("nan")
                bed_score = 0

            tx = transcripts.get(tid)
            blocks = []
            if tx is not None and t_stop_1 >= t_start_1 > 0:
                blocks = transcript_to_genomic_blocks(tx, t_start_1 - 1, t_stop_1)
                if not chrom:
                    chrom = tx.chrom
                if strand not in ("+", "-"):
                    strand = tx.strand
            if not blocks:
                g_start = row.get("ORF_gstart") or row.get("ORF_genomic_start") or ""
                g_stop = row.get("ORF_gstop") or row.get("ORF_genomic_stop") or ""
                try:
                    g_start_i = int(g_start)
                    g_stop_i = int(g_stop)
                    if g_stop_i < g_start_i:
                        g_start_i, g_stop_i = g_stop_i, g_start_i
                    blocks = [(g_start_i - 1, g_stop_i)]
                except (TypeError, ValueError):
                    continue

            orf_id_raw = row.get("ORF_ID") or f"ribocode|{tid}|{t_start_1}-{t_stop_1}"
            rows.append(
                {
                    "orf_id": f"ribocode|{orf_id_raw}",
                    "transcript_id": tid,
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "aa_length": aa_len,
                    "orf_type": orf_type,
                    "raw_score": pval if not math.isnan(pval) else "",
                    "bed_score": bed_score,
                    "blocks": blocks,
                }
            )
    return rows, resolved


_RIBOTISH_GENPOS_RE = re.compile(r"^([^:]+):(\\d+)-(\\d+):([+\\-])\$")


def _parse_ribotish_genpos(s):
    if not s:
        return None
    m = _RIBOTISH_GENPOS_RE.match(s.strip())
    if not m:
        return None
    # GenomePos is 0-based half-open (BED-style).
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def parse_ribotish(path, transcripts, fields):
    resolved = {"score": "", "orf_type": "", "aa_length": ""}
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows, resolved
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            tid = row.get("Tid", "") or ""
            gid = row.get("Gid", "") or (transcripts[tid].gene_id if tid in transcripts else "")

            aa_raw = pick(row, fields, resolved, "aa_length")
            try:
                aa_len = int(aa_raw) if aa_raw else 0
            except (TypeError, ValueError):
                aa_len = 0

            orf_type = pick(row, fields, resolved, "orf_type") or ""

            score_raw = pick(row, fields, resolved, "score")
            try:
                pval = float(score_raw) if score_raw else float("nan")
            except (TypeError, ValueError):
                pval = float("nan")
            bed_score = clamp1000((1.0 - pval) * 1000) if not math.isnan(pval) else 0

            gp = _parse_ribotish_genpos(row.get("GenomePos", ""))
            if gp is None:
                continue
            chrom, start, end, strand = gp

            # ribotish predict reports one genomic span (GenomePos), not
            # per-exon blocks, so the exon structure is recovered by
            # intersecting that span with the transcript's exons.
            tx = transcripts.get(tid)
            if tx is None:
                blocks = [(start, end)]
            else:
                blocks = [(max(start, gs), min(end, ge)) for gs, ge in tx.exons if min(end, ge) > max(start, gs)]
                if not blocks:
                    blocks = [(start, end)]

            rows.append(
                {
                    "orf_id": f"ribotish|{tid}|{chrom}:{blocks[0][0]}-{blocks[-1][1]}:{strand}",
                    "transcript_id": tid,
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "aa_length": aa_len,
                    "orf_type": orf_type,
                    "raw_score": pval if not math.isnan(pval) else "",
                    "bed_score": bed_score,
                    "blocks": blocks,
                }
            )
    return rows, resolved


def _ribotricer_span_from_id(orf_id):
    if not orf_id:
        return None
    parts = orf_id.rsplit("_", 3)
    if len(parts) != 4:
        return None
    try:
        gstart = int(parts[1])
        gend = int(parts[2])
    except ValueError:
        return None
    if gend < gstart:
        gstart, gend = gend, gstart
    return (gstart - 1, gend)


def _ribotricer_blocks_from_id(orf_id, transcripts):
    span = _ribotricer_span_from_id(orf_id)
    if span is None:
        return []
    gstart_0, gend_excl = span
    tid = orf_id.rsplit("_", 3)[0] if orf_id.count("_") >= 3 else ""
    tx = transcripts.get(tid)
    if tx is not None and tx.exons:
        blocks = []
        for ex_gs, ex_ge in tx.exons:
            lo = max(gstart_0, ex_gs)
            hi = min(gend_excl, ex_ge)
            if hi > lo:
                blocks.append((lo, hi))
        if blocks:
            return blocks
    return [(gstart_0, gend_excl)]


def parse_ribotricer(path, transcripts, fields):
    resolved = {"score": "", "orf_type": "", "length": ""}
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows, resolved
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return rows, resolved
        for row in reader:
            status = (row.get("status") or "").lower()
            if status and status != "translating":
                continue
            orf_id_raw = row.get("ORF_ID") or ""
            tid = row.get("transcript_id") or ""
            gid = row.get("gene_id") or ""
            chrom = row.get("chrom") or ""
            strand = row.get("strand") or "+"

            aa_len = aa_from_length(row, fields, resolved)
            orf_type = pick(row, fields, resolved, "orf_type") or ""

            score_raw = pick(row, fields, resolved, "score")
            try:
                pscore = float(score_raw) if score_raw else 0.0
                bed_score = clamp1000(pscore * 1000)
            except (TypeError, ValueError):
                bed_score = 0

            blocks = parse_intervals(row.get("coordinate") or "", seps=("-",))
            if not blocks:
                blocks = _ribotricer_blocks_from_id(orf_id_raw, transcripts)
            if not blocks:
                continue

            rows.append(
                {
                    "orf_id": f"ribotricer|{orf_id_raw}",
                    "transcript_id": tid,
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "aa_length": aa_len,
                    "orf_type": orf_type,
                    "raw_score": score_raw if score_raw else "",
                    "bed_score": bed_score,
                    "blocks": blocks,
                }
            )
    return rows, resolved


def parse_rpbp(path, transcripts, fields):
    """Rp-Bp's predicted-orfs BED extends BED12 with metric columns. The
    file ships with a `#`-prefixed header line; we strip comment lines and
    bind RPBP_COLUMNS as the fieldnames for csv.DictReader so users can
    reference columns by name (e.g. --score-field bayes_factor_var).
    """
    resolved = {"score": "", "length": ""}
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows, resolved
    with open_text(path) as fh:
        lines = (line for line in fh if line and not line.startswith("#"))
        reader = csv.DictReader(lines, fieldnames=RPBP_COLUMNS, delimiter="\\t")
        for row in reader:
            chrom = row.get("seqname", "")
            try:
                start = int(row.get("start") or 0)
            except ValueError:
                continue
            name = row.get("id", "")
            strand = row.get("strand", "+") if row.get("strand") in ("+", "-") else "+"
            try:
                block_count = int(row.get("num_exons") or 0)
                block_sizes = [int(x) for x in (row.get("exon_lengths") or "").rstrip(",").split(",") if x]
                block_starts = [
                    int(x) for x in (row.get("exon_genomic_relative_starts") or "").rstrip(",").split(",") if x
                ]
            except ValueError:
                continue
            if not block_sizes or len(block_sizes) != block_count:
                continue
            blocks = [(start + bs, start + bs + sz) for bs, sz in zip(block_starts, block_sizes)]
            blocks.sort()

            aa_len = aa_from_length(row, fields, resolved)

            score_raw = pick(row, fields, resolved, "score")
            try:
                score_val = float(score_raw) if score_raw else float("nan")
            except (TypeError, ValueError):
                score_val = float("nan")
            bed_score = clamp1000(min(score_val, 30.0) * 1000.0 / 30.0) if not math.isnan(score_val) else 0

            orf_type = "canonical"
            tid = name.rsplit("_", 1)[0] if "_" in name else name
            gid = transcripts[tid].gene_id if tid in transcripts else ""

            rows.append(
                {
                    "orf_id": f"rpbp|{name}",
                    "transcript_id": tid,
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "aa_length": aa_len,
                    "orf_type": orf_type,
                    "raw_score": score_val if not math.isnan(score_val) else "",
                    "bed_score": bed_score,
                    "blocks": blocks,
                }
            )
    return rows, resolved


_PRICE_LOCATION_RE = re.compile(r"^(?P<chrom>.+?)(?P<strand>[+-]):(?P<blocks>.+)\$")


def _parse_price_location(loc):
    if not loc:
        return "", "+", []
    m = _PRICE_LOCATION_RE.match(loc.strip())
    if not m:
        return "", "+", []
    blocks = []
    for tok in m.group("blocks").split("|"):
        tok = tok.strip()
        if not tok or "-" not in tok:
            continue
        a, b = tok.split("-", 1)
        try:
            start = int(a)
            end = int(b)
        except ValueError:
            continue
        if end < start:
            start, end = end, start
        blocks.append((start, end))
    blocks.sort()
    return m.group("chrom"), m.group("strand"), blocks


def _price_score_from_p(p):
    if p is None or p <= 0:
        return 1000
    if p >= 1:
        return 0
    return clamp1000(min(-math.log10(p) * 100, 1000))


def _price_transcript_id(orf_id_raw, orf_type):
    """Recover the transcript id from PRICE's `Id` column.

    PRICE builds `Id` as `<transcript_id>_<orf_type>_<n>` (an increasing
    counter within the transcript). The split anchors on the row's own
    `orf_type` value so transcript ids that themselves contain underscores
    (true for some non-Ensembl annotations) survive intact.
    """
    if orf_type:
        m = re.match(r"^(.+)_" + re.escape(orf_type) + r"_\\d+\$", orf_id_raw)
        if m:
            return m.group(1)
    parts = orf_id_raw.rsplit("_", 2)
    if len(parts) == 3 and parts[2].isdigit():
        return parts[0]
    return orf_id_raw.split("_", 1)[0]


def parse_price(path, transcripts, fields):
    resolved = {"score": "", "orf_type": ""}
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows, resolved
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return rows, resolved
        for row in reader:
            orf_id_raw = (row.get("Id") or "").strip()
            if not orf_id_raw:
                continue
            chrom, strand, blocks = _parse_price_location(row.get("Location") or "")
            if not blocks:
                continue

            orf_type = (pick(row, fields, resolved, "orf_type") or "").strip()

            tid = _price_transcript_id(orf_id_raw, orf_type)
            t = transcripts.get(tid)
            if t is not None:
                gid = t.gene_id
            else:
                # PRICE's `Gene` column concatenates every gene overlapping the
                # ORF's genomic span with an underscore, so this fallback value
                # may itself be a multi-gene string; passed through as-is,
                # downstream joins against a single-gene matrix won't match it.
                gid = (row.get("Gene") or "").strip()

            length_nt = sum(b[1] - b[0] for b in blocks)
            aa_len = max(0, (length_nt - 3) // 3) if length_nt > 0 else 0

            score_raw = pick(row, fields, resolved, "score")
            try:
                pval = float(score_raw) if score_raw else None
            except (TypeError, ValueError):
                pval = None
            bed_score = _price_score_from_p(pval)

            rows.append(
                {
                    "orf_id": f"price|{orf_id_raw}",
                    "transcript_id": tid,
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "aa_length": aa_len,
                    "orf_type": orf_type,
                    "raw_score": pval if pval is not None else "",
                    "bed_score": bed_score,
                    "blocks": blocks,
                }
            )
    return rows, resolved


PARSERS = {
    "ribocode": parse_ribocode,
    "ribotish": parse_ribotish,
    "ribotricer": parse_ribotricer,
    "rpbp": parse_rpbp,
    "price": parse_price,
}


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------


def write_versions():
    with open("versions.yml", "w") as fh:
        yaml.safe_dump(
            {"${task.process}": {"python": platform.python_version()}},
            fh,
            default_flow_style=False,
            sort_keys=False,
        )


def main():
    if CALLER not in PARSERS:
        sys.exit(f"orfnormalise: unknown caller='{CALLER}'; expected one of {sorted(PARSERS)}")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-field", default=None)
    parser.add_argument("--orf-type-field", default=None)
    parser.add_argument("--length-field", default=None)
    parser.add_argument("--aa-length-field", default=None)
    parser.add_argument(
        "--smorf-max-aa",
        type=int,
        default=100,
        help="Maximum aa_length flagged as a small ORF in `is_smorf` (default: 100)",
    )
    args = parser.parse_args(shlex.split("${args}"))
    if args.smorf_max_aa < 1:
        sys.exit(f"orfnormalise: --smorf-max-aa must be >= 1, got {args.smorf_max_aa}")

    fields = {
        "score": _resolve_chain(CALLER, "score", args.score_field),
        "orf_type": _resolve_chain(CALLER, "orf_type", args.orf_type_field),
        "length": _resolve_chain(CALLER, "length", args.length_field),
        "aa_length": _resolve_chain(CALLER, "aa_length", args.aa_length_field),
    }

    transcripts = load_transcripts(GTF)
    rows, resolved_columns = PARSERS[CALLER](INPUT, transcripts, fields)

    bed_lines = []
    tsv_rows = []
    seen = set()
    unmapped = defaultdict(int)

    for r in rows:
        orf_id = r["orf_id"]
        if orf_id in seen:
            continue
        seen.add(orf_id)
        orf_type_native = (r["orf_type"] or "").strip()
        orf_class, matched = classify(CALLER, orf_type_native)
        if not matched:
            unmapped[orf_type_native] += 1
        aa_length = int(r["aa_length"])
        blocks = r["blocks"]
        bed_lines.append(emit_bed12(r["chrom"], blocks, orf_id, int(r["bed_score"]), r["strand"]))
        tsv_rows.append(
            emit_tsv_row(
                orf_id=orf_id,
                caller=CALLER,
                sample_id=SAMPLE_ID,
                chrom=r["chrom"],
                start=blocks[0][0],
                end=blocks[-1][1],
                strand=r["strand"],
                gene_id=r["gene_id"],
                transcript_id=r["transcript_id"],
                orf_class=orf_class,
                aa_length=aa_length,
                score=r["raw_score"],
                orf_type_native=orf_type_native,
                is_smorf=0 < aa_length <= args.smorf_max_aa,
            )
        )

    # An unmapped token would fall through to `other` and pick the wrong merge strategy.
    if unmapped:
        details = ", ".join(f"{token!r} ({count})" for token, count in sorted(unmapped.items()))
        sys.exit(f"orfnormalise: unmapped ORF type value(s) for caller '{CALLER}': {details}")

    write_outputs(OUT_BED, OUT_TSV, bed_lines, tsv_rows, resolved_columns, unmapped, args.smorf_max_aa)
    write_versions()


main()
