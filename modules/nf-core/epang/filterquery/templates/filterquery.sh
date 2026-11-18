#!/usr/bin/env bash
set -euo pipefail

# EPA-ng drops the reference alignment columns that are all gaps, and aborts the whole run with
# "does not appear to have any non-gap sites" if a query has no residue in the columns left.
# A column is "occupied" when at least one reference sequence has a residue there. A query is kept
# when it has at least min_query_columns residues in occupied columns, those residues cover at
# least min_reference_coverage of the occupied columns, and at least one of them is unambiguous:
# EPA-ng also reads '?', N, X and O in either case as gaps, and aborts on a query with nothing else
# in the columns it keeps.
#
# What counts as a residue depends on the alphabet, which is taken from the reference: nucleotide
# when at least 90% of its non-gap characters are A, C, G, T, U or N, otherwise protein. '-' and '.'
# are gaps in both. In nucleotide sequences N and '?' are gaps too, since they carry no information,
# while in protein sequences N is asparagine and counts as a residue.
#
# Both alignments must have the same width. The query file is only written when a query is kept,
# since an empty query file makes EPA-ng fail.

export LC_ALL=C

zcat -f reference/* > reference.fasta
zcat -f query/* > query.fasta

nucleotide=\$(awk '
    /^>/ { next }
    { sub(/\\r\$/, ""); x = \$0; all += gsub(/[^-.?]/, "", x); x = \$0; nuc += gsub(/[ACGTUNacgtun]/, "", x) }
    END { print ((all > 0 && nuc / all >= 0.9) ? 1 : 0) }
' reference.fasta)
echo "EPANG_FILTERQUERY: \$([ "\$nucleotide" = 1 ] && echo nucleotide || echo protein) alphabet"

awk -v nucleotide="\$nucleotide" -v min_cols=${min_query_columns} -v min_cov=${min_reference_coverage} -v table=${prefix}.filterquery.tsv -v kept_fasta=kept.fasta '
    function fail(msg) { print "EPANG_FILTERQUERY: " msg > "/dev/stderr"; failed = 1; exit 1 }
    function flush_ref(   nd, nd2, i, d, d2) {
        if (!in_rec) { return }
        if (width == 0) { width = length(seq); for (j = 1; j <= width; j++) { todo[j] = 1; todo2[j] = 1 } }
        else if (length(seq) != width) { fail("reference sequence " name " has " length(seq) " columns, expected " width) }
        nd = 0
        for (j in todo) { if (substr(seq, j, 1) !~ gapn) { d[++nd] = j } }
        for (i = 1; i <= nd; i++) { delete todo[d[i]] }
        nd2 = 0
        for (j in todo2) { if (substr(seq, j, 1) !~ /[-.?NnXxOo]/) { d2[++nd2] = j } }
        for (i = 1; i <= nd2; i++) { delete todo2[d2[i]] }
        in_rec = 0
    }
    function flush_query(   k, n, m, cov, keep) {
        if (!in_rec) { return }
        if (seq != "" && length(seq) != width) { fail("query " name " has " length(seq) " columns, the reference has " width) }
        n = 0; m = 0
        if (seq != "") {
            for (k = 1; k <= ncols; k++) { if (substr(seq, cols[k], 1) !~ gapn) { n++ } }
            for (k = 1; k <= ncols2; k++) { if (substr(seq, cols2[k], 1) !~ /[-.?NnXxOo]/) { m++ } }
        }
        cov = (ncols > 0) ? n / ncols : 0
        keep = (n >= min_cols && cov >= min_cov && m >= 1)
        printf "%s\\t%d\\t%d\\t%.6f\\t%d\\t%s\\n", name, ncols, n, cov, m, (keep ? "true" : "false") > table
        if (keep) { printf ">%s\\n%s\\n", header, seq > kept_fasta }
        in_rec = 0
    }
    BEGIN { gapn = (nucleotide == 1) ? "[-.?Nn]" : "[-.]"; printf "query\\treference_columns\\tquery_columns\\treference_coverage\\tunambiguous_columns\\tkept\\n" > table }
    FILENAME != current {
        current = FILENAME
        if (current == "query.fasta") {
            flush_ref()
            for (j = 1; j <= width; j++) {
                if (!(j in todo)) { cols[++ncols] = j }
                if (!(j in todo2)) { cols2[++ncols2] = j }
            }
        }
    }
    { sub(/\\r\$/, "") }
    /^>/ {
        if (current == "reference.fasta") { flush_ref() } else { flush_query() }
        header = substr(\$0, 2); name = header; sub(/[ \\t].*/, "", name); seq = ""; in_rec = 1
        next
    }
    { seq = seq \$0 }
    END { if (!failed && current == "query.fasta") { flush_query() } }
' reference.fasta query.fasta

if [ -s kept.fasta ]; then
    gzip -nc kept.fasta > ${prefix}.query.fasta.gz
fi
rm -f reference.fasta query.fasta kept.fasta
