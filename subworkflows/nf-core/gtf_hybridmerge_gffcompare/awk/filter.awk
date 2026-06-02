# Class-code subsetter for a gffcompare *.annotated.gtf
#
# Reads a GTF whose transcript rows carry `class_code "X";` (set by
# gffcompare against the reference annotation), keeps only the records
# belonging to transcripts whose class_code is in the caller-supplied
# `codes_csv` (passed in via gawk's -v from the subworkflow), and emits
# the survivors in their original input order.
#
# Buffer-then-emit so the order of transcript vs exon rows in the input
# does not matter: keep_tx[] is fully populated before we print anything.

BEGIN {
    n = split(codes_csv, arr, /,[[:space:]]*/)
    for (i = 1; i <= n; i++) keep_code[arr[i]] = 1
}

# Drop comments and dot-strand rows (gffcompare emits the latter for
# ambiguous transcripts that have no useful class assignment).
/^#/      { next }
$7 == "." { next }

# Buffer every record that carries a transcript_id; record which tids
# pass the class-code filter (only transcript rows carry class_code).
{
    if (match($0, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if (match($0, /class_code "([^"]+)"/, c) && c[1] in keep_code) {
            keep_tx[tid] = 1
        }
        lines[NR]    = $0
        line_tid[NR] = tid
    }
}

END {
    # Iterate the buffer in numeric-index order so the output preserves
    # input order; gawk's PROCINFO knob gives us that ordering.
    PROCINFO["sorted_in"] = "@ind_num_asc"
    for (i in lines) {
        if (line_tid[i] in keep_tx) print lines[i]
    }
}
