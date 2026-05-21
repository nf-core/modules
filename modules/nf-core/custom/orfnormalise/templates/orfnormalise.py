#!/usr/bin/env python3
"""Normalise an ORF caller's per-sample output to a unified BED12 + sidecar TSV.

One template, five callers. \${caller} (from meta.caller) selects the
per-caller parser, classifier and score mapper.

Harmonised orf_class vocabulary written into the sidecar TSV:
    canonical_cds: ORF maps to an annotated CDS (including truncated /
                   extended variants of one).
    uORF:          upstream ORF (5'UTR-resident).
    dORF:          downstream ORF (3'UTR-resident).
    novel_u:       novel / intergenic ORF not assigned to an annotated CDS.
    smORF:         small ORF (aa_length <= 100). Promoted regardless of
                   location-based class so downstream tools can treat smORFs
                   uniformly.
    other:         internal / overlap / frame variants and anything else.

SCORE_DIRECTIONS records whether a caller's native confidence score is
lower-is-better ('min', e.g. a p-value) or higher-is-better ('max', e.g.
a Bayes factor or phase score). The merger uses this to pick the best
per-ORF call; this module uses it implicitly when mapping the raw score
to the BED 0-1000 score column.
"""
import csv
import gzip
import io
import math
import platform
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

csv.field_size_limit(sys.maxsize)

CALLER = "${caller}"
INPUT = Path("${orfs_table}")
GTF = Path("${gtf}")
SAMPLE_ID = "${sample_id}"
OUT_BED = Path("${prefix}.normalised.bed12")
OUT_TSV = Path("${prefix}.normalised.tsv")

TSV_HEADER = (
    "orf_id\\tcaller\\tsample_id\\tchrom\\tstart\\tend\\tstrand\\t"
    "gene_id\\ttranscript_id\\torf_class\\taa_length\\tscore"
)

SCORE_DIRECTIONS = {
    "ribocode": "min",
    "ribotish": "min",
    "ribotricer": "max",
    "rpbp": "max",
    "price": "min",
}

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


def open_text(path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
    return open(p, "r", encoding="utf-8")


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
    return "\\t".join([
        chrom, str(start), str(end), name, str(score), strand,
        str(start), str(end), "0",
        str(block_count), block_sizes, block_starts,
    ])


def emit_tsv_row(orf_id, caller, sample_id, chrom, start, end, strand,
                 gene_id, transcript_id, orf_class, aa_length, score):
    return "\\t".join([
        orf_id, caller, sample_id, chrom, str(start), str(end), strand,
        gene_id, transcript_id, orf_class, str(aa_length), str(score),
    ])


def write_outputs(bed_path, tsv_path, bed_lines, tsv_rows):
    bed_lines = [l for l in bed_lines if l]
    with open(bed_path, "w") as bh:
        for l in bed_lines:
            bh.write(l + "\\n")
    with open(tsv_path, "w") as th:
        th.write(TSV_HEADER + "\\n")
        for r in tsv_rows:
            th.write(r + "\\n")


def reclassify_smorf(orf_class, aa_length):
    if isinstance(aa_length, int) and 0 < aa_length <= 100:
        return "smORF"
    return orf_class


# ----------------------------------------------------------------------------
# Per-caller classifiers
# ----------------------------------------------------------------------------

def classify_ribotish(tis_type):
    if not tis_type:
        return "other"
    t = tis_type.lower()
    if "5'utr" in t or "uorf" in t:
        return "uORF"
    if "3'utr" in t or "dorf" in t:
        return "dORF"
    if "annotated" in t and "truncated" not in t and "extended" not in t:
        return "canonical_cds"
    if "extended" in t or "truncated" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_ribocode(orf_type):
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "uorf" in t or "5'utr" in t:
        return "uORF"
    if "dorf" in t or "3'utr" in t:
        return "dORF"
    if "annotated" in t or "ccds" in t:
        return "canonical_cds"
    if "internal" in t:
        return "other"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_ribotricer(orf_type):
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "uorf" in t:
        return "uORF"
    if "dorf" in t:
        return "dORF"
    if "annotated" in t or "ccds" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_rpbp(orf_type):
    if not orf_type:
        return "other"
    t = orf_type.lower()
    if "five_prime" in t or "uorf" in t:
        return "uORF"
    if "three_prime" in t or "dorf" in t:
        return "dORF"
    if "canonical" in t or "annotated" in t:
        return "canonical_cds"
    if "novel" in t or "intergenic" in t:
        return "novel_u"
    return "other"


def classify_price(orf_type):
    if not orf_type:
        return "other"
    t = orf_type.strip()
    if t in ("CDS", "Ext", "Trunc", "Variant"):
        return "canonical_cds"
    if t in ("uORF", "uoORF"):
        return "uORF"
    if t == "dORF":
        return "dORF"
    if t == "ncRNA":
        return "novel_u"
    return "other"


CLASSIFIERS = {
    "ribocode": classify_ribocode,
    "ribotish": classify_ribotish,
    "ribotricer": classify_ribotricer,
    "rpbp": classify_rpbp,
    "price": classify_price,
}


# ----------------------------------------------------------------------------
# Per-caller parsers
#
# Each parser returns a DataFrame with columns:
#   orf_id, transcript_id, gene_id, chrom, strand, aa_length, orf_type,
#   raw_score, bed_score, blocks
# `blocks` is a list of (genomic_start_0, genomic_end_excl) tuples in
# genomic order. `bed_score` is the integer 0-1000 BED score. `raw_score`
# is the native caller score (p-value, phase score, Bayes factor) kept for
# the sidecar TSV.
# ----------------------------------------------------------------------------

EMPTY_PARSED = pd.DataFrame(
    columns=["orf_id", "transcript_id", "gene_id", "chrom", "strand",
             "aa_length", "orf_type", "raw_score", "bed_score", "blocks"]
)


def parse_ribocode(path, transcripts):
    if not path.exists() or path.stat().st_size == 0:
        return EMPTY_PARSED.copy()
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return EMPTY_PARSED.copy()
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
            try:
                orf_length_nt = int(row.get("ORF_length") or row.get("ORFlength") or 0)
            except ValueError:
                orf_length_nt = 0
            # RiboCode's predicted_orfs.txt has ORF_length (nt) + AAseq but no
            # AA_length column; derive aa_length nominally as (nt - 3) / 3.
            aa_len = max(0, (orf_length_nt - 3) // 3) if orf_length_nt > 0 else 0
            orf_type = row.get("ORF_type") or row.get("Type") or ""
            try:
                pval = float(row.get("pval_combined") or row.get("Pval_combined") or "1")
                bed_score = max(0, min(1000, int(round((1.0 - pval) * 1000))))
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
            rows.append({
                "orf_id": f"ribocode|{orf_id_raw}",
                "transcript_id": tid,
                "gene_id": gid,
                "chrom": chrom,
                "strand": strand,
                "aa_length": aa_len,
                "orf_type": orf_type,
                "raw_score": pval if pval == pval else "",
                "bed_score": bed_score,
                "blocks": blocks,
            })
    return pd.DataFrame(rows) if rows else EMPTY_PARSED.copy()


_RIBOTISH_GENPOS_RE = re.compile(r"^([^:]+):(\\d+)-(\\d+):([+\\-])\$")


def _parse_ribotish_genpos(s):
    if not s:
        return None
    m = _RIBOTISH_GENPOS_RE.match(s.strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)) - 1, int(m.group(3)), m.group(4)


def _parse_ribotish_blocks(s):
    if not s or s in ("-", "."):
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            a, b = tok.split(":", 1)
        elif "-" in tok:
            a, b = tok.split("-", 1)
        else:
            continue
        try:
            a_i, b_i = int(a), int(b)
        except ValueError:
            continue
        if b_i < a_i:
            a_i, b_i = b_i, a_i
        out.append((a_i - 1, b_i))
    out.sort()
    return out


def parse_ribotish(path, transcripts):
    if not path.exists() or path.stat().st_size == 0:
        return EMPTY_PARSED.copy()
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            tid = row.get("Tid", "") or ""
            gid = row.get("Gid", "") or (transcripts[tid].gene_id if tid in transcripts else "")
            try:
                aa_len = int(row.get("AALen", "0") or 0)
            except ValueError:
                aa_len = 0
            orf_type = row.get("TisType", "") or row.get("TISType", "")
            # Ribo-TISH predict emits TISPvalue, RiboPvalue, FisherPvalue (and
            # their Qvalue equivalents). Fields are written as the literal
            # string "None" when ribotish couldn't compute a value, so we
            # have to fall back to the next available column rather than
            # bailing on a parse error. Prefer FisherPvalue (combined TIS
            # + frame), then RiboPvalue, then TISPvalue.
            pval = float("nan")
            for col in ("FisherPvalue", "Pvalcombined", "RiboPvalue", "TISPvalue", "Pvalue"):
                raw = row.get(col)
                if raw is None or raw == "" or raw == "None":
                    continue
                try:
                    pval = float(raw)
                    break
                except (TypeError, ValueError):
                    continue
            bed_score = max(0, min(1000, int(round((1.0 - pval) * 1000)))) if pval == pval else 0

            gp = _parse_ribotish_genpos(row.get("GenomePos", ""))
            if gp is None:
                continue
            chrom, start, end, strand = gp

            blocks = _parse_ribotish_blocks(row.get("Blocks", ""))
            if not blocks:
                tx = transcripts.get(tid)
                if tx is None:
                    blocks = [(start, end)]
                else:
                    blocks = []
                    for gs, ge in tx.exons:
                        lo = max(start, gs)
                        hi = min(end, ge)
                        if hi > lo:
                            blocks.append((lo, hi))
                    if not blocks:
                        blocks = [(start, end)]

            rows.append({
                "orf_id": f"ribotish|{tid}|{chrom}:{blocks[0][0]}-{blocks[-1][1]}:{strand}",
                "transcript_id": tid,
                "gene_id": gid,
                "chrom": chrom,
                "strand": strand,
                "aa_length": aa_len,
                "orf_type": orf_type,
                "raw_score": pval if pval == pval else "",
                "bed_score": bed_score,
                "blocks": blocks,
            })
    return pd.DataFrame(rows) if rows else EMPTY_PARSED.copy()


def _parse_ribotricer_coord(s):
    if not s:
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok or "-" not in tok:
            continue
        a, b = tok.split("-", 1)
        try:
            a_i, b_i = int(a), int(b)
        except ValueError:
            continue
        if b_i < a_i:
            a_i, b_i = b_i, a_i
        out.append((a_i - 1, b_i))
    out.sort()
    return out


def _ribotricer_span_from_id(orf_id):
    """Parse the genomic span out of a ribotricer ORF_ID of shape
    '<transcript_id>_<gstart>_<gend>_<length_nt>'. Returns (gstart_0,
    gend_excl) or None on parse failure. ORF_ID coordinates are
    1-based-inclusive per ribotricer convention; we return 0-based
    half-open to match BED12.
    """
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
    """Derive ORF block structure from a ribotricer ORF_ID.

    ribotricer's `detect-orfs` output does not carry the intron-aware
    `coordinate` column from `prepare-orfs` - the ORF_ID encodes only the
    outer genomic span. To recover proper multi-exon blocks we intersect
    that span with the host transcript's exon structure from the GTF,
    matching the fallback ribotish uses. Falls back to a single block
    when no GTF / transcript lookup is available.
    """
    span = _ribotricer_span_from_id(orf_id)
    if span is None:
        return []
    gstart_0, gend_excl = span
    # Try to recover the transcript_id from the ORF_ID prefix.
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


def parse_ribotricer(path, transcripts):
    if not path.exists() or path.stat().st_size == 0:
        return EMPTY_PARSED.copy()
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return EMPTY_PARSED.copy()
        for row in reader:
            status = (row.get("status") or "").lower()
            if status and status != "translating":
                continue
            orf_id_raw = row.get("ORF_ID") or ""
            tid = row.get("transcript_id") or ""
            gid = row.get("gene_id") or ""
            chrom = row.get("chrom") or ""
            strand = row.get("strand") or "+"
            try:
                length_nt = int(row.get("length") or 0)
            except ValueError:
                length_nt = 0
            aa_len = max(0, (length_nt - 3) // 3) if length_nt > 0 else 0
            orf_type = row.get("ORF_type") or ""
            try:
                pscore = float(row.get("phase_score") or "0")
                bed_score = max(0, min(1000, int(round(pscore * 1000))))
            except (TypeError, ValueError):
                bed_score = 0

            blocks = _parse_ribotricer_coord(row.get("coordinate") or "")
            if not blocks:
                blocks = _ribotricer_blocks_from_id(orf_id_raw, transcripts)
            if not blocks:
                continue

            rows.append({
                "orf_id": f"ribotricer|{orf_id_raw}",
                "transcript_id": tid,
                "gene_id": gid,
                "chrom": chrom,
                "strand": strand,
                "aa_length": aa_len,
                "orf_type": orf_type,
                "raw_score": "",
                "bed_score": bed_score,
                "blocks": blocks,
            })
    return pd.DataFrame(rows) if rows else EMPTY_PARSED.copy()


def parse_rpbp(path, transcripts):
    if not path.exists() or path.stat().st_size == 0:
        return EMPTY_PARSED.copy()
    rows = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\\t")
            if len(cols) < 12:
                continue
            chrom = cols[0]
            try:
                start = int(cols[1])
            except ValueError:
                continue
            name = cols[3]
            strand = cols[5] if cols[5] in ("+", "-") else "+"
            try:
                block_count = int(cols[9])
                block_sizes = [int(x) for x in cols[10].rstrip(",").split(",") if x]
                block_starts = [int(x) for x in cols[11].rstrip(",").split(",") if x]
            except ValueError:
                continue
            if not block_sizes or len(block_sizes) != block_count:
                continue
            blocks = [(start + bs, start + bs + sz) for bs, sz in zip(block_starts, block_sizes)]
            blocks.sort()

            # Rp-Bp's predicted-orfs BED extends the BED12 with metric
            # columns (orf_num, orf_len, p_translated_mean/var,
            # p_background_mean/var, bayes_factor_mean/var, chi_square_p,
            # x_{1,2,3}_sum, profile_sum). No orf_type column - the
            # final-prediction-set BED is the curated post-filter output, so
            # all rows are confident translated ORFs; default the class to
            # canonical_cds and let GTF lookup (downstream merger) refine.
            try:
                orf_len_nt = int(cols[13]) if len(cols) > 13 else 0
            except ValueError:
                orf_len_nt = 0
            try:
                bf_mean = float(cols[18]) if len(cols) > 18 else float("nan")
            except ValueError:
                bf_mean = float("nan")
            orf_type = "canonical"
            aa_len = max(0, (orf_len_nt - 3) // 3) if orf_len_nt > 0 else 0

            if bf_mean == bf_mean:
                bed_score = max(0, min(1000, int(round(min(bf_mean, 30.0) * 1000.0 / 30.0))))
            else:
                bed_score = 0

            tid = name.rsplit("_", 1)[0] if "_" in name else name
            gid = transcripts[tid].gene_id if tid in transcripts else ""

            rows.append({
                "orf_id": f"rpbp|{name}",
                "transcript_id": tid,
                "gene_id": gid,
                "chrom": chrom,
                "strand": strand,
                "aa_length": aa_len,
                "orf_type": orf_type,
                "raw_score": bf_mean if bf_mean == bf_mean else "",
                "bed_score": bed_score,
                "blocks": blocks,
            })
    return pd.DataFrame(rows) if rows else EMPTY_PARSED.copy()


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
    s = int(round(min(-math.log10(p) * 100, 1000)))
    return max(0, min(1000, s))


def parse_price(path, transcripts):
    if not path.exists() or path.stat().st_size == 0:
        return EMPTY_PARSED.copy()
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        if reader.fieldnames is None:
            return EMPTY_PARSED.copy()
        for row in reader:
            orf_id_raw = (row.get("Id") or "").strip()
            if not orf_id_raw:
                continue
            chrom, strand, blocks = _parse_price_location(row.get("Location") or "")
            if not blocks:
                continue
            orf_type = (row.get("Type") or "").strip()
            # PRICE Id format is `<transcript_id>_<type>_<index>`; take the
            # leading underscore-delimited token as the host transcript.
            tid = orf_id_raw.split("_", 1)[0]
            gid = (row.get("Gene") or "").strip()
            if not gid:
                t = transcripts.get(tid)
                if t is not None:
                    gid = t.gene_id

            length_nt = sum(b[1] - b[0] for b in blocks)
            aa_len = max(0, (length_nt - 3) // 3) if length_nt > 0 else 0

            try:
                pval = float(row.get("p value") or row.get("p_value") or "")
            except (TypeError, ValueError):
                pval = None
            bed_score = _price_score_from_p(pval)

            rows.append({
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
            })
    return pd.DataFrame(rows) if rows else EMPTY_PARSED.copy()


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


def main():
    if CALLER not in PARSERS:
        sys.exit(f"orfnormalise: unknown meta.caller='{CALLER}'; expected one of {sorted(PARSERS)}")

    transcripts = load_transcripts(GTF)
    df = PARSERS[CALLER](INPUT, transcripts)

    bed_lines = []
    tsv_rows = []
    seen = set()
    classify = CLASSIFIERS[CALLER]

    for _, r in df.iterrows():
        orf_id = r["orf_id"]
        if orf_id in seen:
            continue
        seen.add(orf_id)
        orf_class = reclassify_smorf(classify(r["orf_type"]), int(r["aa_length"]))
        blocks = r["blocks"]
        bed_lines.append(
            emit_bed12(r["chrom"], blocks, orf_id, int(r["bed_score"]), r["strand"])
        )
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
                aa_length=int(r["aa_length"]),
                score=r["raw_score"],
            )
        )

    write_outputs(OUT_BED, OUT_TSV, bed_lines, tsv_rows)
    write_versions()


main()
