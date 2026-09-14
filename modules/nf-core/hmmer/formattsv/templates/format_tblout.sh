#!/usr/bin/env bash
set -euo pipefail

# Turns HMMER's --tblout format into a clean TSV, one row per hit. Every column
# but the last (a free-text description) is whitespace-separated. The real
# column count (N) is read per file from its dash/alignment header row: every
# token there is a run of dashes, so counting them is unambiguous, unlike the
# column-name row, whose labels (e.g. "target name") also contain spaces and
# don't map one-to-one onto columns. Each data row's first N-1 whitespace
# tokens become the named columns; the rest is rejoined into the final column.
#
# Whether HMMER leaves a space between "#" and the first column's dashes
# depends on that column's width. With no gap, awk's default splitting glues
# them into one token ("#---------"), so the first column needs separate
# handling below.

labels=(${labels.join(' ')})
files=(${files.join(' ')})

{
    printf 'profile\ttarget_name\ttarget_accession\tquery_name\tquery_accession\tfull_evalue\tfull_score\tfull_bias\tbest_domain_evalue\tbest_domain_score\tbest_domain_bias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription\n'

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
} | gzip -n > ${prefix}.tblout.tsv.gz
