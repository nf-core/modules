process MOTIFMATCHR_MATCHMOTIFS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-motifmatchr:1.32.0--r45ha27e39d_0':
        'quay.io/biocontainers/bioconductor-motifmatchr:1.32.0--r45ha27e39d_0' }"

    input:
    tuple val(meta), path(peaks), path(fasta), path(motifs)

    output:
    tuple val(meta), path("*.motif_matches.mtx"), emit: matches
    tuple val(meta), path("*.motifmatchr_peaks.tsv"), emit: peaks
    tuple val(meta), path("*.motifmatchr_motifs.tsv"), emit: motifs
    tuple val("${task.process}"), val('motifmatchr'), eval("Rscript -e 'cat(as.character(utils::packageVersion(\"motifmatchr\")))'"), topic: versions, emit: versions_motifmatchr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript - $peaks $fasta $motifs $prefix $args <<'RSCRIPT'
    `%||%` <- function(x, y) if (is.null(x)) y else x
    arguments <- commandArgs(trailingOnly = TRUE)
    if (length(arguments) < 4L) stop('Expected peaks, FASTA, motifs, and output prefix')
    peaks_file <- arguments[[1L]]
    fasta_file <- arguments[[2L]]
    motifs_file <- arguments[[3L]]
    out_prefix <- arguments[[4L]]
    extra_arguments <- arguments[-seq_len(4L)]
    p_cutoff_index <- match('--p-cutoff', extra_arguments)
    p_cutoff <- if (is.na(p_cutoff_index)) 5e-5 else as.numeric(extra_arguments[[p_cutoff_index + 1L]])
    if (is.na(p_cutoff) || p_cutoff <= 0 || p_cutoff >= 1) stop('--p-cutoff must be between 0 and 1')

    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(IRanges)
        library(Matrix)
        library(motifmatchr)
        library(Rsamtools)
        library(TFBSTools)
    })

    peaks <- read.delim(peaks_file, header = FALSE, comment.char = '#', stringsAsFactors = FALSE)
    if (ncol(peaks) < 3L || nrow(peaks) == 0L) stop('--peaks must contain at least three BED columns and one region')
    region_id <- if (ncol(peaks) >= 4L) as.character(peaks[[4L]]) else rep('', nrow(peaks))
    missing_id <- is.na(region_id) | !nzchar(region_id)
    region_id[missing_id] <- paste0('peak_', which(missing_id))
    region_id <- make.unique(region_id)
    regions <- GRanges(
        seqnames = peaks[[1L]],
        ranges = IRanges(start = as.integer(peaks[[2L]]) + 1L, end = as.integer(peaks[[3L]])),
        region_id = region_id
    )

    motif_list <- readJASPARMatrix(motifs_file, matrixClass = 'PFM')
    if (length(motif_list) == 0L) stop('No position frequency matrices could be read from --motifs')
    matches <- matchMotifs(motif_list, regions, genome = FaFile(fasta_file), out = 'matches', p.cutoff = p_cutoff)
    match_matrix <- motifMatches(matches)
    colnames(match_matrix) <- make.unique(colnames(match_matrix))

    writeMM(match_matrix, paste0(out_prefix, '.motif_matches.mtx'))
    write.table(
        data.frame(row_index = seq_len(length(regions)), region_id = region_id, chrom = as.character(seqnames(regions)), start = start(regions) - 1L, end = end(regions)),
        paste0(out_prefix, '.motifmatchr_peaks.tsv'), quote = FALSE, row.names = FALSE, sep = '\t'
    )
    write.table(
        data.frame(column_index = seq_len(ncol(match_matrix)), motif_id = colnames(match_matrix)),
        paste0(out_prefix, '.motifmatchr_motifs.tsv'), quote = FALSE, row.names = FALSE, sep = '\t'
    )
    RSCRIPT
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%%MatrixMarket matrix coordinate pattern general\\n1 1 0\\n' > ${prefix}.motif_matches.mtx
    printf 'row_index\\tregion_id\\tchrom\\tstart\\tend\\n1\\tpeak_1\\tchr1\\t0\\t1\\n' > ${prefix}.motifmatchr_peaks.tsv
    printf 'column_index\\tmotif_id\\n1\\tmotif_1\\n' > ${prefix}.motifmatchr_motifs.tsv
    """
}
