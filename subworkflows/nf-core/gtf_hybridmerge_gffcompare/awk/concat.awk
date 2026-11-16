# Backbone + novel concat with synthesised parent gene rows
#
# Reads two files in order: the canonical backbone GTF (full
# gene/transcript/exon structure) followed by the filtered novel GTF
# (gffcompare's annotated output downstream of the class-code filter,
# which contains only transcript/exon/CDS rows - gffcompare strips
# parent gene rows). Emits a single hybrid GTF in which every novel
# gene_id that was absent from the backbone gets a synthesised `gene`
# row whose coordinates span the union of all child transcripts.

BEGIN { FS = "\t"; OFS = "\t" }

# Comments are dropped from both files.
/^#/ { next }

# Backbone gene rows: record the gid as already covered so we never
# synthesise a duplicate later. The row itself goes through the buffer
# as-is.
$3 == "gene" {
    if (match($9, /gene_id "([^"]+)"/, m)) gene_seen[m[1]] = 1
}

# Transcript rows: if this is a novel gene_id (no backbone gene row),
# track the union span and the chromosome/source/strand of the first
# sighting so END can emit a coherent synth gene row.
$3 == "transcript" {
    if (match($9, /gene_id "([^"]+)"/, m)) {
        gid = m[1]
        if (!(gid in gene_seen)) {
            if (!(gid in min_start) || $4 < min_start[gid]) min_start[gid] = $4
            if (!(gid in max_end)   || $5 > max_end[gid])   max_end[gid]   = $5
            if (!(gid in tx_chr)) {
                tx_chr[gid]    = $1
                tx_src[gid]    = $2
                tx_strand[gid] = $7
            }
        }
    }
}

# Every record gets buffered. Synth gene rows are inserted in END at the
# point where the first transcript for each novel gene_id appears.
{ lines[++nlines] = $0 }

END {
    PROCINFO["sorted_in"] = "@ind_num_asc"
    for (i in lines) {
        line = lines[i]
        n = split(line, f, "\t")
        if (n >= 8 && f[3] == "transcript" && match(line, /gene_id "([^"]+)"/, m)) {
            gid = m[1]
            if (!(gid in gene_seen) && !(gid in synth)) {
                # Synth gene row uses the union span across all novel
                # transcripts of this gid - matters for tools (e.g.
                # genome browsers, span-based filters) that take gene
                # coords at face value rather than re-deriving from
                # children.
                print tx_chr[gid], tx_src[gid], "gene", min_start[gid], max_end[gid], ".", tx_strand[gid], ".", "gene_id \"" gid "\";"
                synth[gid] = 1
            }
        }
        print line
    }
}
