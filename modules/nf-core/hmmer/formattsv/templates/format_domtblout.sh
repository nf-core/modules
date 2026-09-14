#!/usr/bin/env bash
set -euo pipefail

# Turns HMMER's --domtblout format into a clean TSV, one row per domain hit.
# Same parsing rule as format_tblout.sh, including the first-column dash-row
# handling: the real column count (N) is read per file from its dash/alignment
# header row, then each data row's first N-1 whitespace tokens become the
# named columns and the rest is rejoined into the final ("description of
# target") column.

labels=(${labels.join(' ')})
files=(${files.join(' ')})

{
    printf 'profile\ttarget_name\ttarget_accession\ttarget_length\tquery_name\tquery_accession\tquery_length\tfull_evalue\tfull_score\tfull_bias\tdomain_number\tdomain_count\tc_evalue\ti_evalue\tdomain_score\tdomain_bias\thmm_from\thmm_to\tali_from\tali_to\tenv_from\tenv_to\tacc\tdescription\n'

    for i in "\${!files[@]}"; do
        label="\${labels[i]}"
        file="\${files[i]}"
        zcat -f "\$file" | awk -v label="\$label" '
            /^#/ {
                if (N == 0) {
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
                    if (ok && n > 0) N = n
                }
                next
            }
            N == 0 { next }
            {
                printf "%s", label
                for (i = 1; i < N; i++) printf "\\t%s", \$i
                rest = \$N
                for (i = N + 1; i <= NF; i++) rest = rest " " \$i
                printf "\\t%s\\n", rest
            }
        '
    done
} | gzip -n > ${prefix}.domtblout.tsv.gz
