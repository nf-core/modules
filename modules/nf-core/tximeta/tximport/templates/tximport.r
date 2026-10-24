#!/usr/bin/env Rscript --vanilla

# Script for importing and processing transcript-level quantifications.
# Written by Lorena Pantano, later modified by Jonathan Manning and Rob Syme,
# and released under the MIT license.
#
# Samples are imported in blocks so that memory scales with the gene-level
# (gene x sample) matrices rather than with several transcript x sample
# matrices. Each block is imported with tximport() and summed to gene level
# with the same operations as tximport::summarizeToGene(), and its
# transcript-level tables are written to part files that are joined at the end.
# Gene-level values that depend on the whole cohort (missing-length filling and
# countsFromAbundance scaling) are computed with tximport's own helpers once all
# blocks are in, so outputs match a single tximport() + summarizeToGene() run.

# tximport::summarizeToGene() reorders gene-level output rows via base R's
# rowsum(reorder=TRUE), which sorts by the process locale's collation order.
# Fix LC_COLLATE so row order (and therefore output content) is reproducible
# across execution environments regardless of their default locale.
Sys.setlocale("LC_COLLATE", "C")

# Loading required libraries
library(tximport)

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(y) scan(text=y, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(y) y[2]), names = lapply(args_vals, function(y) y[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Read Transcript Metadata from a Given Path
#'
#' This function reads transcript metadata from a specified file path. The file is expected to
#' be a tab-separated values file with a header, containing transcript information. Columns are
#' selected by name using the tx_col, gene_id_col, and gene_name_col parameters. The function
#' checks if the file is empty and stops execution with an error message if so. Additional
#' processing is done to match the quantified transcripts, including adding missing entries and
#' reordering based on the quantified transcript IDs.
#'
#' @param tinfo_path The file path to the transcript information file.
#' @param tx_ids Transcript IDs in the order they appear in the quantification files.
#' @param tx_col Column name for transcript IDs.
#' @param gene_id_col Column name for gene IDs.
#' @param gene_name_col Column name for gene names.
#'
#' @return A list containing five elements:
#' - `transcript`: A data frame with transcript IDs, gene IDs, and gene names, indexed by transcript IDs.
#' - `gene`: A data frame with unique gene IDs and gene names.
#' - `tx2gene`: A data frame mapping transcript IDs to gene IDs.
#' - `tx2gene_augmented`: The full tx2gene table actually used by tximport (input
#'   mappings plus self-mappings appended for any transcripts present in the
#'   quantification output but missing from the input tx2gene), with the input
#'   column headers preserved.
#' - `extra`: A character vector of transcript IDs found in quantification output but missing from the tx2gene file.

read_transcript_info <- function(tinfo_path, tx_ids, tx_col, gene_id_col, gene_name_col){
    info <- file.info(tinfo_path)
    if (info\$size == 0) {
        stop("tx2gene file is empty")
    }

    # Read file with actual header to handle variable column counts correctly
    raw_info <- read.csv(tinfo_path, sep="\t", header = TRUE, check.names = FALSE)

    # Resolve column references (use named columns if present, otherwise positional)
    resolved_tx_col        <- if (tx_col        %in% colnames(raw_info)) tx_col        else colnames(raw_info)[1]
    resolved_gene_id_col   <- if (gene_id_col   %in% colnames(raw_info)) gene_id_col   else colnames(raw_info)[2]
    resolved_gene_name_col <- if (gene_name_col %in% colnames(raw_info)) gene_name_col else colnames(raw_info)[3]

    transcript_info <- data.frame(
        tx        = raw_info[[resolved_tx_col]],
        gene_id   = raw_info[[resolved_gene_id_col]],
        gene_name = raw_info[[resolved_gene_name_col]],
        check.names = FALSE
    )

    extra <- setdiff(tx_ids, as.character(transcript_info[["tx"]]))
    if (length(extra) > 0) {
        warning(
            length(extra), " transcripts found in quantification output but missing from ",
            "the tx2gene mapping (GTF). These will be included in transcript-level outputs ",
            "but excluded from gene-level summaries. This usually means the transcript FASTA ",
            "and GTF are from different sources or versions. First 5: ",
            paste(head(extra, 5), collapse = ", ")
        )
    }
    transcript_info <- rbind(transcript_info, data.frame(tx=extra, gene_id=extra, gene_name=extra, check.names = FALSE))
    transcript_info <- transcript_info[match(tx_ids, transcript_info[["tx"]]), ]
    rownames(transcript_info) <- transcript_info[["tx"]]

    # Restore input column headers so the augmented file can be fed back to tximport
    tx2gene_augmented <- transcript_info
    colnames(tx2gene_augmented) <- c(resolved_tx_col, resolved_gene_id_col, resolved_gene_name_col)

    list(transcript = transcript_info,
        gene = unique(transcript_info[,2:3]),
        tx2gene = transcript_info[,1:2],
        tx2gene_augmented = tx2gene_augmented,
        extra = extra)
}

#' Import transcript-level quantifications for a set of samples
#'
#' @param files Named vector of quantification files (names are sample names).
#' @param quant_type Quantification type: 'salmon', 'kallisto' or 'rsem'.
#'
#' @return The tximport() result, with transcript x sample matrices `abundance`, `counts` and `length`.

import_quants <- function(files, quant_type) {
    if (quant_type == "rsem") {
        tximport(files, type = 'rsem', txIn = TRUE, txOut = TRUE)
    } else {
        # Inferential replicates never reach the outputs, so don't read them
        tximport(files, type = quant_type, txOut = TRUE, dropInfReps = TRUE)
    }
}

#' Split long-double row sums into columns of doubles
#'
#' rowSums() and rowMeans() accumulate in long double, so a running total kept
#' in double precision between blocks would round differently from a single
#' rowMeans() over all samples. This returns, for each row, a few doubles whose
#' sum is exactly the row's long-double total. Passing them back through
#' rowSums() ahead of the next block's columns continues the accumulation from
#' exactly where it stopped.
#'
#' @param x A numeric matrix.
#'
#' @return A matrix with one row per row of x and as many columns as needed.

rowsum_parts <- function(x) {
    parts <- x[, 0, drop = FALSE]
    repeat {
        residual <- rowSums(cbind(x, -parts))
        if (!any(residual != 0, na.rm = TRUE)) {
            return(parts)
        }
        parts <- cbind(parts, residual)
    }
}

#' Row means over n samples from the parts returned by rowsum_parts()
#'
#' Pads the parts with zero columns up to n so that rowMeans() performs the
#' same long-double division by n as it would on the full matrix. Rows are
#' processed in chunks to keep the padding small.
#'
#' @param parts Matrix returned by rowsum_parts() over all n samples.
#' @param n Number of samples.
#'
#' @return Numeric vector identical to rowMeans() of the full matrix.

rowmeans_from_parts <- function(parts, n) {
    means <- numeric(nrow(parts))
    chunk <- max(1L, 1e6 %/% n)
    for (start in seq(1L, nrow(parts), by = chunk)) {
        rows <- start:min(nrow(parts), start + chunk - 1L)
        padding <- matrix(0, length(rows), n - ncol(parts))
        means[rows] <- rowMeans(cbind(parts[rows, , drop = FALSE], padding))
    }
    means
}

#' Write identifier columns and a value matrix as a tab-separated table
#'
#' Rows are written in chunks so that no full-size data frame copy of the
#' matrix is made. The output is identical to a single write.table() call.
#'
#' @param file Output file path.
#' @param values Matrix with sample names as column names.
#' @param ids Data frame of identifier columns to put first, one row per
#'   written row of values, or NULL to write values only.
#' @param rows Rows of values to write.

write_tsv <- function(file, values, ids = NULL, rows = seq_len(nrow(values))) {
    con <- file(file, "w")
    on.exit(close(con))
    chunk <- 20000L
    for (start in seq(1L, max(length(rows), 1L), by = chunk)) {
        idx <- seq.int(start, length.out = min(chunk, length(rows) - start + 1L))
        table <- as.data.frame(values[rows[idx], , drop = FALSE], optional = TRUE)
        if (!is.null(ids)) {
            table <- data.frame(ids[idx, , drop = FALSE], table, check.names = FALSE)
        }
        write.table(table, con, sep="\t", quote=FALSE, row.names = FALSE, col.names = start == 1L)
    }
}

################################################
################################################
## Main script starts here                    ##
################################################
################################################

# Set defaults and classes
opt <- list(
    tx_col = "transcript_id",
    gene_id_col = "gene_id",
    gene_name_col = "gene_name",
    block_size = "100"
)

# Apply parameter overrides from ext.args
args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (ao %in% names(opt)) {
        opt[[ao]] <- args_opt[[ao]]
    }
}
block_size <- suppressWarnings(as.integer(opt\$block_size))
if (is.na(block_size) || block_size < 1) {
    stop("--block_size must be a positive integer, got: ", opt\$block_size)
}

# Define pattern for file names based on quantification type
if ('$quant_type' == "rsem") {
    # RSEM: .isoforms.results files are flat in quants/ directory (not in subdirectories)
    fns <- list.files('quants', pattern = 'isoforms\\\\.results\$', recursive = TRUE, full.names = TRUE)
    names(fns) <- gsub("\\\\.isoforms\\\\.results\$", "", basename(fns))
} else {
    # Salmon/Kallisto: files are in sample subdirectories
    pattern <- ifelse('$quant_type' == "kallisto", "abundance.tsv", "quant.sf")
    fns <- list.files('quants', pattern = pattern, recursive = TRUE, full.names = TRUE)
    names(fns) <- basename(dirname(fns))
}
if (length(fns) == 0) {
    stop("No quantification files found for quant_type '$quant_type'")
}
if (anyDuplicated(names(fns))) {
    stop("Duplicate sample names: ", paste(unique(names(fns)[duplicated(names(fns))]), collapse = ", "))
}

prefix <- ''
if ('$task.ext.prefix' != 'null'){
    prefix = '$task.ext.prefix'
} else if ('$meta.id' != 'null'){
    prefix = '$meta.id'
}

# Transcript-level tables are written block by block into part files, which
# are joined column-wise at the end
transcript_tables <- c(
    abundance = "transcript_tpm.tsv",
    counts = "transcript_counts.tsv",
    length = "transcript_lengths.tsv"
)
part_dir <- "tximport_parts"
dir.create(part_dir, showWarnings = FALSE)
part_files <- list()
part_types <- list()

blocks <- split(seq_along(fns), ceiling(seq_along(fns) / block_size))
for (b in seq_along(blocks)) {
    samples <- blocks[[b]]
    txi <- import_quants(fns[samples], '$quant_type')
    message("Imported block ", b, " of ", length(blocks), " (", length(samples), " samples)")

    if (b == 1) {
        tx_ids <- rownames(txi[["abundance"]])

        # Read transcript and gene metadata
        transcript_info <- read_transcript_info('$tx2gene', tx_ids, opt\$tx_col, opt\$gene_id_col, opt\$gene_name_col)

        # Gene grouping as tximport::summarizeToGene() builds it from a tx2gene
        # table with one row per quantified transcript
        gene_id <- factor(transcript_info\$tx2gene[["gene_id"]])

        # Gene x sample accumulators. Like tximport's own matrices, these start
        # as logical NA so that their type follows the data (integer or double).
        new_gene_matrix <- function() {
            matrix(NA, nlevels(gene_id), length(fns), dimnames = list(levels(gene_id), names(fns)))
        }
        gene_abundance <- new_gene_matrix()
        gene_counts <- new_gene_matrix()
        gene_weighted_length <- new_gene_matrix()
        length_sum_parts <- NULL
    } else if (!identical(rownames(txi[["abundance"]]), tx_ids)) {
        stop("Transcript IDs in block ", b, " differ from those in the first block")
    }

    # Gene-level sums as in tximport::summarizeToGene(). Each sample's column
    # is summed independently, so a block gives the same values as the cohort.
    gene_abundance[, samples] <- rowsum(txi[["abundance"]], gene_id)
    gene_counts[, samples] <- rowsum(txi[["counts"]], gene_id)
    gene_weighted_length[, samples] <- rowsum(txi[["abundance"]] * txi[["length"]], gene_id)

    # Running long-double transcript length totals, for rowMeans() at the end
    length_sum_parts <- rowsum_parts(cbind(length_sum_parts, txi[["length"]]))

    for (assay in names(transcript_tables)) {
        part_file <- file.path(part_dir, sprintf("%s.%05d.tsv", assay, b))
        ids <- if (b == 1) transcript_info\$transcript[, 1:2] else NULL
        write_tsv(part_file, txi[[assay]], ids)
        part_files[[assay]] <- c(part_files[[assay]], part_file)
        part_types[[assay]] <- c(part_types[[assay]], typeof(txi[[assay]]))
    }
    rm(txi)
    invisible(gc())
}

# Join the transcript-level part files into the final tables
for (assay in names(transcript_tables)) {
    files <- part_files[[assay]]

    # tximport stores a matrix that holds only whole numbers as integer, which
    # prints differently from double (100000 rather than 1e+05). If other blocks
    # made the cohort-wide matrix double, rewrite such blocks as double.
    if ("double" %in% part_types[[assay]]) {
        for (i in which(part_types[[assay]] != "double")) {
            part <- read.delim(files[i], check.names = FALSE, colClasses = "character", quote = "")
            id_cols <- if (i == 1) 1:2 else integer(0)
            values <- as.matrix(part[, setdiff(seq_along(part), id_cols), drop = FALSE])
            storage.mode(values) <- "double"
            write_tsv(files[i], values, if (i == 1) part[, id_cols] else NULL)
        }
    }

    output_file <- paste0(prefix, ".", transcript_tables[[assay]])
    if (length(files) == 1) {
        file.rename(files, output_file)
    } else if (system2("paste", shQuote(files), stdout = output_file) != 0) {
        stop("Failed to join transcript-level part files into ", output_file)
    }
}
unlink(part_dir, recursive = TRUE)

# Finish the gene-level summary exactly as tximport::summarizeToGene() does
gene_length <- gene_weighted_length / gene_abundance
rm(gene_weighted_length)
ave_length_samp <- rowmeans_from_parts(length_sum_parts, length(fns))
ave_length_samp_gene <- tapply(ave_length_samp, gene_id, mean)
stopifnot(all(names(ave_length_samp_gene) == rownames(gene_length)))
gene_length <- tximport:::replaceMissingLength(gene_length, ave_length_samp_gene)

# Remove fake gene entries created by unmapped transcripts (where gene_id
# was set to the transcript ID). These would break downstream processes
# like SummarizedExperiment that try to match gene IDs against the tx2gene.
real_genes <- which(!rownames(gene_abundance) %in% transcript_info\$extra)

gene_info <- transcript_info\$gene[match(rownames(gene_abundance)[real_genes], transcript_info\$gene[["gene_id"]]),]
rownames(gene_info) <- NULL

write_gene_table <- function(values, suffix) {
    write_tsv(paste0(prefix, ".", suffix), values, gene_info, real_genes)
}
write_gene_table(gene_length, "gene_lengths.tsv")
write_gene_table(gene_abundance, "gene_tpm.tsv")
write_gene_table(gene_counts, "gene_counts.tsv")

# Counts from abundance are scaled with cohort-wide totals (and, for
# lengthScaledTPM, each gene's mean length across all samples), using all
# genes including those removed above, as summarizeToGene() does
scaled_tables <- c(lengthScaledTPM = "gene_counts_length_scaled.tsv", scaledTPM = "gene_counts_scaled.tsv")
for (cfa in names(scaled_tables)) {
    scaled_counts <- tximport:::makeCountsFromAbundance(gene_counts, gene_abundance, gene_length, cfa)
    write_gene_table(scaled_counts, scaled_tables[[cfa]])
    rm(scaled_counts)
}

write.table(transcript_info\$tx2gene_augmented,
    paste0(prefix, ".tx2gene_augmented.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(prefix, "R_sessionInfo.log", sep = '.'))
citation("tximeta")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
tximeta.version <- as.character(packageVersion('tximeta'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-tximeta:', tximeta.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
