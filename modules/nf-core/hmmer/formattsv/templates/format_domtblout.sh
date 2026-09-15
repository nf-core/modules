#!/usr/bin/env bash
set -euo pipefail

# Turns HMMER's --domtblout format into a clean TSV, one row per domain hit.
# Same parsing rule as format_tblout.sh, including the first-column dash-row
# handling: the real column count (N) is read per file from its dash/alignment
# header row, then each data row's first N-1 whitespace tokens become the
# named columns and the rest is rejoined into the final ("description of
# target") column.
#
# Labels/files are quoted per element rather than bare-joined: an unquoted
# join word-splits any label or path containing whitespace and misaligns the
# label/file pairing for every entry after it.
#
# Once a header is found, a data row is only ever treated as a trailing
# comment (hmmsearch's "# Program:"/"# Date:"/"# [ok]" run-metadata footer)
# once the exact sentinel line "#" (nothing else) has been seen -- a target
# name that happens to start with "#" is otherwise valid data, not a comment.

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
            N > 0 && \$0 == "#" { intrailer = 1 }
            intrailer { next }
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
