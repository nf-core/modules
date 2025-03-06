#!/usr/bin/env Rscript

# Calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
# Author: Ira Iosub

suppressPackageStartupMessages(library(riboWaltz))
suppressPackageStartupMessages(library(dplyr))

# =========
# Options and paths
# =========

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

# =========
# Functions for riboWaltz analysis
# =========

#' Export P-site offset table produced by the `riboWaltz::psite` function for each BAM file based on sample name
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df dataframe produced by the riboWaltz psite function with P-site offsets for each read-length and a sample column
#'
#' @return filtered dataframe for the sample specified with sample_name. It creates and saves a TSV file to disk.

export_offsets <- function(sample_name, df) {

    df <- dplyr::filter(df, sample == sample_name)
    data.table::fwrite(df, paste0(getwd(), "/", sample_name, ".offset.tsv.gz"), sep = "\t")
    return(df)

}

#' Export P-sites table produced by the `riboWaltz::psite_info` function for each BAM file based on sample name
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df_list dataframe list produced by the riboWaltz psite_info function with P-sites information for each alignment
#'
#' @return filtered dataframe for the sample specified with sample_name. It creates and saves a TSV file to disk.

export_psites <- function(sample_name, df_list) {

    df <- df_list[[sample_name]]
    df <- dplyr::mutate(df, sample = sample_name)
    data.table::fwrite(df, paste0(getwd(), "/", sample_name, ".psite.tsv.gz"), sep = "\t")
    return(data.table::data.table(df))

}

#' Export codon P-site and RPF coverage tables produced by the `riboWaltz::codon_coverage` function for each BAM file based on sample name
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df.ls dataframe list produced by the riboWaltz psite_info function with P-sites information for all samples
#' @param annotation.df annotation dataframe list produced by the riboWaltz create_annotation function with CDS, UTR length etc info for each transcript
#'
#' @return This function does not return a value. It creates and saves a TSV file to disk.

export_codon_coverage_tables <- function(sample_name, df.ls, annotation.df) {

    rpf_coverage.dt <- riboWaltz::codon_coverage(df.ls, annotation = annotation.df, sample = sample_name, psite = FALSE)
    psite_coverage.dt <- riboWaltz::codon_coverage(df.ls, annotation = annotation.df, sample = sample_name, psite = TRUE)
    data.table::fwrite(rpf_coverage.dt, paste0(getwd(),"/",sample_name, ".codon_coverage_rpf.tsv.gz"), sep = "\t")
    data.table::fwrite(psite_coverage.dt, paste0(getwd(),"/",sample_name, ".codon_coverage_psite.tsv.gz"), sep = "\t")

}

#' Export CDS P-site count tables produced by the `riboWaltz::cds_coverage` function for each BAM file based on sample name, for the entire CDS and a defined CDS window
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param cds_coverage dataframe produced by the riboWaltz cds_coverage function with P-sites counts over the entire CDS of transcripts
#' @param cds_window_coverage dataframe produced by the riboWaltz cds_coverage function with P-sites counts over the defined CDS window of transcripts using the start_nts and stop_nts params
#' @param exclude_start number of nucleotides from start codon that were excluded when defining the CDS window for counting (provide same calue as start_nts param)
#' @param exclude_stop number of nucleotides from stop codon that were excluded when defining the CDS window for counting (provide same calue as stop_nts param)
#'
#' @return This function does not return a value. It creates and saves two TSV files to disk.

export_cds_coverage_tables <- function(sample_name, cds_coverage, cds_window_coverage, exclude_start, exclude_stop) {

    cols <- c("transcript", "length_cds", sample_name)
    cds_coverage.dt <- cds_coverage[,..cols]
    cds_window_coverage.dt <- cds_window_coverage[,..cols]
    # Export CDS coverage tables
    data.table::fwrite(cds_coverage.dt, paste0(getwd(),"/",sample_name, ".cds_coverage_psite.tsv.gz"), sep = "\t")
    data.table::fwrite(cds_window_coverage.dt, paste0(getwd(), "/", sample_name, ".cds_","plus", exclude_start, "nt_minus", exclude_stop, "nt_coverage_psite.tsv.gz"), sep = "\t")

}

#' Export plots read-length distribution of reads used for P-site offset identification produced by the `riboWaltz::rlength_distr` function
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df_list dataframe list containig data for all samples
#'
#' @return This function does not return a value. It creates and saves a PDF file to disk.

plot_length_bins <- function(sample_name, df_list) {

    comparison_list <- list()
    comparison_list[["start_codon"]] <- df_list[[sample_name]][end5 <= cds_start & end3 >= cds_start]
    comparison_list[["whole_sample"]] <- df_list[[sample_name]]

    if(nrow(comparison_list[["start_codon"]]) == 0) {
        comparison_list <- list()
        comparison_list[["whole_sample"]] <- df_list[[sample_name]]
        rpf_list <- list("All" = c("whole_sample"))
        length_dist_split <-  riboWaltz::rlength_distr(comparison_list,
                                            sample = rpf_list,
                                            multisamples = "average",
                                            plot_style = "split",
                                            colour = c("gray70"))

    } else {

        rpf_list <- list("Only_start" = c("start_codon"), "All" = c("whole_sample"))

        length_dist_split <-  riboWaltz::rlength_distr(comparison_list,
                                            sample = rpf_list,
                                            multisamples = "average",
                                            plot_style = "facet",
                                            colour = c("#699FC4", "gray70"))
    }

    length_dist_split.gg <- length_dist_split[["plot"]] +
        ggplot2::ggtitle(sample_name)


    ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".length_bins_for_psite.pdf"), length_dist_split.gg, dpi = 400, width = 10, height = 5)

}

#' Export meta-heatmaps of read extremities around start and stop codons produced by the `riboWaltz::rends_heat` function
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df_list dataframe list produced by the riboWaltz bam_to_list function
#' @param annotation annotation dataframe list produced by the riboWaltz create_annotation function with CDS, UTR length etc info for each transcript
#'
#' @return This function does not return a value. It creates and saves a PDF file to disk.

plot_metaheatmap <- function(sample_name, df_list, annotation) {

    ends_heatmap <- riboWaltz::rends_heat(df_list, annotation, sample = sample_name,
                                cl = 100, utr5l = 25, cdsl = 40, utr3l = 25)

    ends_heatmap.gg <- ends_heatmap[[paste0("plot_", sample_name)]] +
        ggplot2::ylim(20, 45)

    ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".ends_heatmap.pdf"), ends_heatmap.gg, dpi = 400, width = 12, height = 8)

}

#' Export meta-heatmaps of read extremities around start and stop codons produced by the `riboWaltz::codon_usage_psite` function
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param psite_info_ls dataframe list produced by the riboWaltz psite_info function with P-sites information for each alignment
#' @param frequency_normalization ogical value whether to normalize the
#' 64 codon usage indexes for the corresponding codon frequencies in coding
#' sequences. Default is TRUE.
#'
#' @return This function does not return a value. It creates and saves a PDF file to disk.

plot_codon_usage <- function(sample_name, psite_info_ls, frequency_normalization = TRUE) {

    psite.ls <- psite_info_ls[sample_name]

    cu_barplot <- riboWaltz::codon_usage_psite(psite.ls, annotation = annotation.dt, sample = sample_name,
                                    fasta_genome = TRUE,
                                    fastapath = opt\$fasta,
                                    gtfpath = opt\$gtf,
                                    frequency_normalization = frequency_normalization)

    cu_barplot.gg <-cu_barplot[[paste0("plot_", sample_name)]]

    ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".codon_usage.pdf"), cu_barplot.gg, dpi = 400, width = 10, height = 7)
}


#' Exclude a sample from the dataframe list if no reads overlap start codon
#'
#' @param sample_name string specifying the sample name associated with the BAM file
#' @param df_list dataframe list produced by the riboWaltz bam_to_list function
#'
#' @return vector with sample names to be excluded

exclude_samples <- function(sample_name, df_list) {

    sample_list <- list()
    exclude <- c()
    sample_list[["start_codon"]] <- df_list[[sample_name]][end5 <= cds_start & end3 >= cds_start]

    if(nrow(sample_list[["start_codon"]]) == 0) {
        message("No reads overlapping start codon. Removing sample from analysis.")

        exclude <- sample_name

    } else {

        message("This sample will not be excluded.")
    }

    return(exclude)
}

#' Quietly Stops the R Session
#'
#' This function stops the R session without printing any error messages to the console.
#' It temporarily suppresses error messages during the stopping process, ensuring that
#' the session halts smoothly without displaying potentially confusing or distracting output.
#'
#' @return This function does not return a value. It simply stops the R session.

stop_quietly <- function() {

    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
}

#' Save Length Distribution Plot
#'
#' This function generates a plot of read length distributions using the `riboWaltz::rlength_distr` function across a given sample and saves it to disk.
#'
#' @param sample_name A string specifying the name of the sample for which to create the plot.
#' @param dt.ls A list of dataframes containing the data for all samples.
#'
#' @details
#' The function first generates a length distribution plot using the `riboWaltz::rlength_distr` function.
#' The plot is saved to a PDF file in the "ribowaltz_qc" directory.
#'
#' @return This function does not return a value. It creates and saves a PDF file to disk.

save_length_distribution_plot <- function(sample_name, dt.ls) {

    # Read lengths averaged across samples
    length_dist <- riboWaltz::rlength_distr(reads.ls, sample = sample_name, multisamples = "average", cl = 99, colour = "grey70")

    length_dist.gg <- length_dist[[paste0("plot_", sample_name)]] +
        ggplot2::theme(legend.position = "none", legend.title=ggplot2::element_blank())

    ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".length_distribution.pdf"), length_dist.gg, dpi = 400)

}

#' Save P-site Region Plot
#'
#' This function generates and saves a plot of P-site region distributions using the `riboWaltz::region_psite` across a given sample,
#' displaying the distribution of reads mapping onto annotated regions.
#'
#' @param sample_name A string specifying the name of the sample for which to create the plot.
#' @param dt.ls A list of dataframes containing data for all samples.
#' @param annotation.df A dataframe containing transcript annotations, which include lengths of
#' the 5' UTR, CDS, and 3' UTR for each transcript.
#'
#' @details
#' This function calculates and visualizes the distribution of P-sites across different transcript
#' regions: 5' UTR, CDS, and 3' UTR. The plot is saved as a PDF in the "ribowaltz_qc" directory.
#'
#' @return This function does not return a value. It generates and saves a PDF file to disk.

save_psite_region_plot <- function(sample_name, dt.ls, annotation.df) {

    psite_region <- riboWaltz::region_psite(dt.ls, annotation = annotation.df, sample = sample_name)
    psite_region.gg <- psite_region[["plot"]] +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

    ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".psite_region.pdf"), psite_region.gg, dpi = 400, width = 10)

}

#' Save Frame Plots
#'
#' This function generates and saves plots displaying the distribution of P-sites across
#' different reading frames and transcript regions using the `riboWaltz::frame_psite_length` and `riboWaltz::frame_psite` functions.
#'
#' @param sample_name A string specifying the name of the sample for which to create the plots.
#' @param dt.ls A data.frame list containing data for all samples.
#' @param annotation.df A dataframe containing transcript annotations, which include lengths of
#' the 5' UTR, CDS, and 3' UTR for each transcript.
#' @param min_length An integer specifying the minimum read length to consider for plotting.
#' @param max_length An integer specifying the maximum read length to consider for plotting.
#'
#' @details
#' The function first generates a stratified frame plot using the `riboWaltz::frame_psite_length` function.
#' This function visualizes the percentage of P-sites falling into each reading frame for each read length,
#' across all transcript regions. The plot is then enhanced with axis and background customizations,
#' and saved as a PDF in the "ribowaltz_qc" directory.
#'
#' Next, a second frame plot is generated using `riboWaltz::frame_psite` which aggregates
#' the distribution of P-sites across all read lengths for each transcript region.
#' The plot is similarly enhanced and saved.
#'
#' @return This function does not return a value. It generates and saves two PDF files to disk.

save_frame_plots <- function(sample_name, dt.ls, annotation.df, min_length, max_length) {

    frames_stratified <- riboWaltz::frame_psite_length(dt.ls, region = "all", sample = sample_name, length_range = min_length:max_length, annotation = annotation.df)
    frames_stratified.gg <- frames_stratified[[paste0("plot_", sample_name)]] +
        ggplot2::scale_y_continuous(limits = c(min_length - 0.5, max_length + 0.5), breaks = seq(min(min_length + ((min_length) %% 2), max_length), max(min_length + ((min_length) %% 2), max_length),
                                                                                                by = max(2, floor((max_length - min_length) / 7))))

    ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".frames_stratified.pdf"), frames_stratified.gg, dpi = 600, height = 9 , width = 12)


    frames <- riboWaltz::frame_psite(dt.ls, region = "all", length_range = min_length:max_length, sample = sample_name, annotation = annotation.df, colour = "grey70")
    frames.gg <- frames[[paste0("plot_", sample_name)]]

    ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".frames.pdf"), frames.gg, dpi = 600, height = 9 , width = 9)

}


#' Save Metaprofile P-site Plot
#'
#' This function generates and saves a metaprofile plot showing the distribution of P-sites
#' around start and stop codons for a specific sample using the `riboWaltz::metaprofile_psite` function.
#'
#' @param sample_name A string specifying the name of the sample for which to create the plot.
#' @param df.ls A list containing data for all samples.
#' @param annotation.df A dataframe containing transcript annotations, including information
#' about UTR and CDS lengths for each transcript.
#'
#' @details
#' This function uses the `riboWaltz::metaprofile_psite` function to compute a metaprofile that
#' aggregates P-site positions across all transcripts for the specified sample. The plot visualizes
#' P-site abundance around start and stop codons, taking into account specified regions of
#' 25 nucleotides before the start, 40 nucleotides in the coding sequence, and 25 nucleotides after the stop.
#'
#' The generated plot saved as a PDF file in the "ribowaltz_qc" directory in a wide format.
#'
#' @return This function does not return a value. It generates and saves a PDF file to disk.

save_metaprofile_psite_plot <- function(sample_name, df.ls, annotation.df) {

    metaprofile <- riboWaltz::metaprofile_psite(df.ls, annotation.df, sample = sample_name,
                                    utr5l = 25, cdsl = 40, utr3l = 25, colour = "black")

    metaprofiles.gg <- metaprofile[[paste0("plot_", sample_name)]]

    ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".metaprofile_psite.pdf"), metaprofiles.gg,
                    dpi = 400, width = 12, height = 6) # save in wide format

}
# =========
# Parse parameters for Nextflow
# =========

# Set defaults and classes
opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    threads = '$task.cpus',
    bam = '$bam',
    gtf = '$gtf',
    fasta = '$fasta',
    length_range = NULL,               # specifies a range of read lengths for P-site identification, formatted as two integers separated by a colon (e.g., '26:34'). If unspecified, all lengths are considered.
    periodicity_threshold = NULL,      # filter out read lengths below specified periodicity threshold for P-site identification. must be a value between 10 and 100. If null, no periodicity filtering is done.
    start = TRUE,                      # use the translation initiation site as reference codon for calculating offsets. If 'FALSE', the second last codon is used instead.
    extremity = 'auto',                # specifies if the offset correction step should be based on 5' extremities ('5end') or 3'extremities ('3end') or 'auto' i.e. the optimal extremity is automatically selected.
    start_nts = 42,                    # number of nt from start codon used to exclude P-sites near initiating ribosome when calculating CDS window P-site counts.
    stop_nts = 27,                     # number of nucleotides from stop codon used to exclude P-sites near terminating ribosome when calculating CDS coverage when calculating CDS window P-site counts.
    frequency_normalization = TRUE     # for codon usage index calculation.
)
opt_types <- lapply(opt, class)

# Parse extra arguments
args_opt <- parse_args('$task.ext.args')

# Apply parameter overrides
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check file inputs are valid
for (file_input in c('bam', 'gtf', 'fasta')){
    if (! is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# Set number of threads for data.table
data.table::setDTthreads(opt\$threads)

# Create folder for plots
dir.create("ribowaltz_qc")

# =========
# Load data and prepare annotation
# =========

# Prepare annotation: for each transcript, obtain total length, 5'UTR, CDS and 3'UTR length, respectively.
annotation.dt <- riboWaltz::create_annotation(opt\$gtf)
data.table::fwrite(annotation.dt,
                    paste0(getwd(), "/", strsplit(opt\$gtf, "gtf")[[1]][1],"transcript_info.tsv.gz"),
                    sep = "\t")

# Load BAM file
bam_name <- opt\$output_prefix
names(bam_name) <- strsplit(basename(opt\$bam), ".bam")[[1]][1]

reads.ls <- riboWaltz::bamtolist(bamfolder = getwd(), annotation = annotation.dt, name_samples = bam_name)

stopifnot(length(reads.ls) == 1) # check that a single bam has been provided and will be analysed

# Order named list alphabetically
reads.ls <- reads.ls[order(names(reads.ls))]


# Get filtered reads: keep only the ones with periodicity evidence based on periodicity_threshold, if provided
if (!is.null(opt\$periodicity_threshold)) {

    filtered.ls <- riboWaltz::length_filter(data = reads.ls,
                                length_filter_mode = "periodicity",
                                periodicity_threshold = opt\$periodicity_threshold)
} else {

    filtered.ls <- reads.ls

}

# Additionally filter them by length, if length_range provided
if (!is.null(opt\$length_range)) {

    min_length <- as.integer(strsplit(opt\$length_range, ":")[[1]][1])
    max_length <- as.integer(strsplit(opt\$length_range, ":")[[1]][2])

    filtered.ls <- riboWaltz::length_filter(data = filtered.ls,
                                length_filter_mode = "custom",
                                length_range = min_length:max_length)
}

# Remove sample if no reads left after filtering steps
filtered.ls <- Filter(function(x) dim(x)[1] > 0, filtered.ls)

# Plot length bins used for P-site assignment using remaining reads
lapply(names(filtered.ls), plot_length_bins, df_list = filtered.ls)

# Filter out sample if no reads pass filtering, and stop analysis if not a single sample passes filering
exclude.ls <- lapply(names(filtered.ls), exclude_samples, df_list = filtered.ls)
filtered.ls <- filtered.ls[!names(filtered.ls) %in% exclude.ls]

if (length(filtered.ls) == 0) {

    message("No sample has reads passing filters for P-site identification. Stopping analysis.")
    stop_quietly()
}

# =========
# Identify P-sites
# =========

message("Calculating P-site offsets and P-site positions defined by the riboWaltz method...")

# Compute P-site offsets: temporary and corrected
psite_offset.dt <- riboWaltz::psite(filtered.ls, flanking = 6, extremity = opt\$extremity, start = opt\$start,
                        txt = TRUE, plot = TRUE, plot_format = "pdf", txt_file = paste0(bam_name, ".best_offset.txt"))

# Save offsets for each sample
lapply(unique(psite_offset.dt\$sample), export_offsets, df = psite_offset.dt)


# Update reads information according to the inferred P-sites
filtered_psite.ls <- riboWaltz::psite_info(filtered.ls, psite_offset.dt, site = "psite",
                                fasta_genome = TRUE, refseq_sep = " ",
                                fastapath = opt\$fasta,
                                gtfpath = opt\$gtf)

sample_name.ls <- unique(names(filtered_psite.ls))

# Save psite info for each sample
lapply(sample_name.ls, export_psites, df_list = filtered_psite.ls)

message("Calculating codon and CDS coverage...")

# 1. codon_coverage called by export_codon_coverage_tables computes the number of read footprints or P-sites mapping on each triplet of annotated coding sequences and UTRs.
lapply(sample_name.ls, export_codon_coverage_tables, df.ls = filtered_psite.ls, annotation.df = annotation.dt)

# 2. Compute the number of P-sites mapping on annotated coding sequences or whole transcripts.
# By default, only in-frame P-sites falling in the full annotated coding sequences are considered.
cds_coverage_psite.dt <- riboWaltz::cds_coverage(filtered_psite.ls, annotation = annotation.dt)

# Additionally compute the number of in-frame P-sites per coding sequences excluding specified number of nucleotides, by default the first 15 codons and the last 10 codons
cds_coverage_psite_window.dt <- riboWaltz::cds_coverage(filtered_psite.ls, annotation = annotation.dt, start_nts = opt\$start_nts, stop_nts = opt\$stop_nts)
lapply(sample_name.ls, export_cds_coverage_tables, cds_coverage = cds_coverage_psite.dt, cds_window_coverage = cds_coverage_psite_window.dt,
            exclude_start = opt\$start_nts, exclude_stop = opt\$stop_nts)

# =========
# Diagnostic plots
# =========

message("Generating diagnostic plots...")

# Define min and max length if not provided as a param
if (is.null(opt\$length_range)) {

    min_rl <- as.integer(min(psite_offset.dt[,"length"]))
    max_rl <- as.integer(max(psite_offset.dt[,"length"]))
}

lapply(names(reads.ls), save_length_distribution_plot, dt.ls = reads.ls)

# Metaheatmaps: the abundance of the 5' and 3' extremity of reads mapping on and around the start and the stop codon of annotated CDSs, stratified by their length.
ends_heatmap.gg.ls <- lapply(names(reads.ls), plot_metaheatmap, df_list = reads.ls, annotation = annotation.dt)

# Ribosome profiling data should highlight the CDS of transcripts as the region with the higher percentage of reads.
lapply(sample_name.ls, save_psite_region_plot, dt.ls = filtered_psite.ls, annotation.df = annotation.dt)

# Compute the percentage of P-sites falling in the three possible translation reading frames for 5’ UTRs, CDSs and 3’ UTRs.
# Plots should show an enrichment of P-sites in the first frame on the coding sequence but not the UTRs, as expected for ribosome protected fragments from protein coding mRNAs.
lapply(sample_name.ls, save_frame_plots, dt.ls = filtered_psite.ls, annotation.df = annotation.dt,
        min_length = min_rl, max_length = max_rl)

# Trinucleotide periodicity along coding sequences: metaprofiles (the merge of single, transcript-specific profiles) based on P-sites mapping around the start and the stop codon of annotated CDSs.
lapply(sample_name.ls, save_metaprofile_psite_plot, df.ls = filtered_psite.ls, annotation.df = annotation.dt)

# Codon usage
lapply(names(filtered_psite.ls), plot_codon_usage, psite_info_ls = filtered_psite.ls, frequency_normalization = opt\$frequency_normalization)

message("riboWaltz analysis succesfully completed!")

# =========
# Export versions
# =========

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
ribowaltz.version <- as.character(packageVersion('riboWaltz'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-ribowaltz:', ribowaltz.version)
    ),
'versions.yml')
