#!/usr/bin/env bash
set -euo pipefail

# Turns HMMER's --domtblout format into a clean TSV, one row per domain hit. Same parsing
# rule as format_tblout.sh: the header row's dashes give the real column count (N), each data
# row's first N-1 tokens become the named columns, and the rest is rejoined into the final
# (description) column.
#
# Labels/files are quoted per element, not bare-joined, so a space in a label or path can't
# misalign labels with files.
#
# A data row seen after the exact "#" sentinel (the start of HMMER's run-metadata footer)
# fails loudly instead of being silently dropped -- it means more than one run's output was
# concatenated into this file.

labels=(${labels.collect { "'" + it.toString().replace("'", "'\\''") + "'" }.join(' ')})
files=(${files.collect { "'" + it.toString().replace("'", "'\\''") + "'" }.join(' ')})

{
    printf 'profile\ttarget_name\ttarget_accession\ttarget_length\tquery_name\tquery_accession\tquery_length\tfull_evalue\tfull_score\tfull_bias\tdomain_number\tdomain_count\tc_evalue\ti_evalue\tdomain_score\tdomain_bias\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\tacc\tdescription\n'

    for i in "\${!files[@]}"; do
        label="\${labels[i]}"
        file="\${files[i]}"
        zcat -f "\$file" | awk -v label="\$label" -v expected=23 '
            { sub(/\\r\$/, "") }
            N == 0 && /^#/ {
                n = 0
                ok = 1
                first = \$1
                sub(/^#/, "", first)
                if (length(first) > 0) {
                    if (first ~ /^-+\$/) { n++ } else { ok = 0 }
                }
                if (ok) {
                    for (i = 2; i <= NF; i++) {
                        if (\$i !~ /^-+\$/) { ok = 0; break }
                        n++
                    }
                }
                if (ok && n > 0) {
                    if (n != expected) {
                        print "hmmer/formattsv: expected " expected " domtblout columns, found " n " -- unsupported layout" > "/dev/stderr"
                        exit 1
                    }
                    N = n
                }
                next
            }
            N > 0 && \$0 == "#" { intrailer = 1; next }
            intrailer {
                if (\$0 != "" && \$0 !~ /^#/) {
                    print "hmmer/formattsv: data row seen after the run-metadata footer -- looks like more than one HMMER run concatenated into one file" > "/dev/stderr"
                    exit 1
                }
                next
            }
            {
                if (N == 0) {
                    print "hmmer/formattsv: data row seen before a column-count header was found" > "/dev/stderr"
                    exit 1
                }
                printf "%s", label
                for (i = 1; i < N; i++) printf "\\t%s", \$i
                rest = \$N
                for (i = N + 1; i <= NF; i++) rest = rest " " \$i
                printf "\\t%s\\n", rest
            }
        '
    done
} | gzip -n > ${prefix}.domtblout.tsv.gz
