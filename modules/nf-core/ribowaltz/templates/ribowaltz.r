#!/usr/bin/env Rscript

# Calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
# Author: Ira Iosub
# Usage: Rscript ribowaltz.r --bams SRX11780887.Aligned.toTranscriptome.out.bam,SRX11780888.Aligned.toTranscriptome.out.bam --gtf Homo_sapiens.GRCh38.111_chr20.gtf --fasta Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz --length_range "26:34" --periodicity_threshold 50 --exclude_start 42 --exclude_stop 27


suppressPackageStartupMessages(library(riboWaltz))

# Allow data.table to use all available threads
data.table::setDTthreads(0)

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


# # Function to print help menu
# print_help_menu <- function() {
#   cat("
# Usage: Rscript ribowaltz.r [options]
# Options:
#   --bams FILE                Path to the BAM file(s). If multiple bams provided, they must be comma-separated. (required)
#   --gtf FILE                 Path to the GTF file. (required)
#   --fasta FILE               Path to the FASTA file. (required)
#   --length_range RANGE       Specify a range of read lengths for P-site identification, formatted as two integers separated by a colon (e.g., '26:34'). If unspecified, all lengths are considered.
#   --periodicity_threshold INT   Filter out read lengths below specified periodicity threshold for P-site identification. Must be a value between 10 and 100. (default: NULL)
#   --exclude_start INT        Number of nucleotides from start codon used to exclude P-sites near initiating ribosome when calculating CDS coverage. (default: 42)
#   --exclude_stop INT         Number of nucleotides from stop codon used to exclude P-sites near terminating ribosome when calculating CDS coverage. (default: 27)
#   -h, --help                 Display this help message and exit.
# ")
#   quit(status = 0)  # Exit after displaying help
# }


# =========
# Parse parameters for Nextflow
# =========

# Set defaults and classes
opt <- list(
    output_prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
    bams = '$bam',
    gtf = '$gtf',
    fasta = '$fasta',
    length_range = NULL,
    periodicity_threshold = NULL,
    exclude_start = 42,               # number of nt from start codon used to exclude P-sites near initiating ribosome when calculating CDS P-site counts.
    exclude_stop = 27
)
opt_types <- lapply(opt, class)


# Parse extra arguments
args_opt <- parse_args('$task.ext.args')

# # Check if help is requested
# if ('-h' %in% args_opt || '--help' %in% args_opt) {
#   print_help_menu()
# }

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
for (file_input in c('bams', 'gtf', 'fasta')){
    if (! is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}


# Print inputs and parsed arguments
cat("BAM file(s): ", opt\$bams, "\n")
cat("GTF file: ", opt\$gtf, "\n")
cat("FASTA file: ", opt\$fasta, "\n")
cat("Length range for filtering: ", opt\$length_range, "\n")
cat("Periodicity threshold: ", opt\$periodicity_threshold, "\n")
cat("Exclude start codons (in nt) up to: ", opt\$exclude_start, "\n")
cat("Exclude stop codons (in nt) up to: ", opt\$exclude_stop, "\n")

# Create folder for plots
dir.create("ribowaltz_qc")


# =========
# Functions for riboWaltz analysis
# =========


export_offsets <- function(name, df) {
  
  df <- df[sample == name]
  
  data.table::fwrite(df, paste0(getwd(), "/", name, ".offset.tsv.gz"), sep = "\t")
  return(df)
  
}

export_psites <- function(name, df_list) {
  
  df <- df_list[[name]]
  df$sample <- name
  
  data.table::fwrite(df, paste0(getwd(), "/", name, ".psite.tsv.gz"), sep = "\t")
  return(df)
  
}

export_codon_coverage_tables <- function(sample_name, df.ls, annotation.df) {
  
  rpf_coverage.dt <- codon_coverage(df.ls, annotation = annotation.df, sample = sample_name, psite = FALSE)
  psite_coverage.dt <- codon_coverage(df.ls, annotation = annotation.df, sample = sample_name, psite = TRUE)
  data.table::fwrite(rpf_coverage.dt, paste0(getwd(),"/",sample_name, ".codon_coverage_rpf.tsv.gz"), sep = "\t")
  data.table::fwrite(psite_coverage.dt, paste0(getwd(),"/",sample_name, ".codon_coverage_psite.tsv.gz"), sep = "\t")
  
}

export_cds_coverage_tables <- function(sample_name, cds_coverage, cds_window_coverage, exclude_start, exclude_stop) {
  
  cols <- c("transcript", "length_cds", sample_name)
  cds_coverage.dt <- cds_coverage[,..cols]
  cds_window_coverage.dt <- cds_window_coverage[,..cols]
  
  # Export CDS coverage tables
  data.table::fwrite(cds_coverage.dt, paste0(getwd(),"/",sample_name, ".cds_coverage_psite.tsv.gz"), sep = "\t")
  data.table::fwrite(cds_window_coverage.dt, paste0(getwd(), "/", sample_name, ".cds_","plus", exclude_start, "nt_minus", exclude_stop, "nt_coverage_psite.tsv.gz"), sep = "\t")
  
}

plot_length_bins <- function(sample_name, df_list) {
  
  comparison_list <- list()
  comparison_list[["start_codon"]] <- df_list[[sample_name]][end5 <= cds_start & end3 >= cds_start]
  comparison_list[["whole_sample"]] <- df_list[[sample_name]]
  
  if(nrow(comparison_list[["start_codon"]]) == 0) {
    
    comparison_list <- list()
    comparison_list[["whole_sample"]] <- df_list[[sample_name]]
    
    rpf_list <- list("All" = c("whole_sample"))
    
    length_dist_split <-  rlength_distr(comparison_list,
                                        sample = rpf_list,
                                        multisamples = "average",
                                        plot_style = "split",
                                        colour = c("gray70"))
    
  } else {
    
    rpf_list <- list("Only_start" = c("start_codon"), "All" = c("whole_sample"))
    
    length_dist_split <-  rlength_distr(comparison_list,
                                        sample = rpf_list,
                                        multisamples = "average",
                                        plot_style = "facet",
                                        colour = c("#699FC4", "gray70"))
  }
  
  length_dist_split.gg <- length_dist_split[["plot"]] +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank()) +
    ggplot2::ggtitle(sample_name)
  
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".length_bins_for_psite.pdf"), length_dist_split.gg, dpi = 400, width = 10, height = 5)
  
}

plot_metaheatmap <- function(sample_name, df_list, annotation) {
  
  ends_heatmap <- rends_heat(df_list, annotation, sample = sample_name, 
                             cl = 100,
                             utr5l = 25, cdsl = 40, utr3l = 25)
  
  ends_heatmap.gg <- ends_heatmap[[paste0("plot_", sample_name)]] +
    ggplot2::ylim(20,45)
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".ends_heatmap.pdf"), ends_heatmap.gg, dpi = 400, width = 12, height = 8)
  
  return(ends_heatmap.gg)
}

plot_codon_usage <- function(sample_name, psite_info_ls) {
  
  psite.ls <- psite_info_ls[sample_name]
  
  cu_barplot <- codon_usage_psite(psite.ls, annotation = annotation.dt, sample = sample_name,
                                  fasta_genome = TRUE, 
                                  fastapath = opt\$fasta,
                                  gtfpath = opt\$gtf,
                                  frequency_normalization = FALSE) 
  
  cu_barplot.gg <-cu_barplot[[paste0("plot_", sample_name)]] +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".codon_usage.pdf"), cu_barplot.gg, dpi = 400, width = 10, height = 7)
}

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

stop_quietly <- function() {
  
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

save_length_distribution_plot <- function(sample_name, dt.ls) {
  
  # Read lengths averaged across samples
  length_dist <- rlength_distr(reads.ls, sample = sample_name, multisamples = "average", cl = 99, colour = "grey70")
  
  length_dist.gg <- length_dist[[paste0("plot_", sample_name)]] +
    ggplot2::theme(legend.position = "none", legend.title=ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = "grey70") +
    ggplot2::scale_color_manual(values = "grey30") +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".length_distribution.pdf"), length_dist.gg, dpi = 400)
  
}

save_psite_region_plot <- function(sample_name, dt.ls, annotation.df) {
  
  psite_region <- region_psite(dt.ls, annotation = annotation.df, sample = sample_name)
  psite_region.gg <- psite_region$plot + 
    ggplot2::theme(plot.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  
  ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".psite_region.pdf"), psite_region.gg, dpi = 400, width = 10)
  
}

# A fundamental characteristic of ribosome profiling data is the trinucleotide periodicity of ribosome footprints along coding sequences. 
# Functions frame_psite_length and frame_psite show if, and to which extent, the identified P-sites results in codon periodicity on the CDS. 
# Both functions compute the percentage of P-sites falling in the three possible translation reading frames for 5’ UTRs, CDSs and 3’ UTRs with one difference: 
# frame_psite_length analyses all read lengths separately and generates a heatmap for each transcript region, while frame_psite processes all reads at once, returning three bar plots.
save_frame_plots <- function(sample_name, dt.ls, annotation.df) {
  
  frames_stratified <- frame_psite_length(dt.ls, region = "all", sample = sample_name, length_range = min_length:max_length, annotation = annotation.df)
  frames_stratified.gg <- frames_stratified[[paste0("plot_", sample_name)]] +
    ggplot2::scale_y_continuous(limits = c(min_length - 0.5, max_length + 0.5), breaks = seq(min(min_length + ((min_length) %% 2), max_length), max(min_length + ((min_length) %% 2), max_length), 
                                                                                             by = max(2, floor((max_length - min_length) / 7))))
  
  ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".frames_stratified.pdf"), frames_stratified.gg, dpi = 600, height = 9 , width = 12)
  
  
  frames <- frame_psite(dt.ls, region = "all", length_range = min_length:max_length, sample = sample_name, annotation = annotation.df, colour = "grey70")
  frames.gg <- frames[[paste0("plot_", sample_name)]] +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/", sample_name, ".frames.pdf"), frames.gg, dpi = 600, height = 9 , width = 9)
  
}

save_metaprofile_psite_plot <- function(sample_name, df.ls, annotation.df) {
  
  metaprofile <- metaprofile_psite(df.ls, annotation.df, sample = sample_name,
                                   utr5l = 25, cdsl = 40, utr3l = 25, colour = "black")
  
  metaprofiles.gg <- metaprofile[[paste0("plot_", sample_name)]]
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".metaprofile_psite.pdf"), metaprofiles.gg,
                  dpi = 400, width = 12, height = 6) # save in wide format
  
}


# =========
# Load data and prepare annotation
# =========

# Prepare annotation: for each transcript, obtain total length, 5'UTR, CDS and 3'UTR length, respectively.
annotation.dt <- create_annotation(opt\$gtf)
data.table::fwrite(annotation.dt, 
                   paste0(getwd(), "/", strsplit(opt\$gtf, "gtf")[[1]][1],"transcript_info.tsv.gz"),
                   sep = "\t")

# Load BAM files
bams <- as.list(strsplit(opt\$bams, ",")[[1]]) # if multiple bams are provided, creates a list

name_of_bams <- lapply(bams, function(x) strsplit(basename(x), ".transcriptome.sorted.bam")[[1]][1])
name_of_bams <- lapply(name_of_bams, function(x) strsplit(basename(x), ".Aligned.toTranscriptome.sorted.out.bam")[[1]][1])
name_of_bams <- lapply(name_of_bams, function(x) strsplit(basename(x), ".Aligned.toTranscriptome.out.bam")[[1]][1])


names(name_of_bams) <- lapply(bams, function(x) strsplit(basename(x), ".bam")[[1]][1])

# Load bams
reads.ls <- bamtolist(bamfolder = getwd(), 
                      annotation = annotation.dt,
                      name_samples = unlist(name_of_bams))

# Order named list alphabetically
reads.ls <- reads.ls[order(names(reads.ls))]


# Get filtered reads: keep only the ones with periodicity evidence based on periodicity_threshold, if privided
if (!is.null(opt\$periodicity_threshold)) {

  filtered.ls <- length_filter(data = reads.ls,
                             length_filter_mode = "periodicity",
                             periodicity_threshold = opt\$periodicity_threshold)
} else {

  filtered.ls <- reads.ls

}

# Additionally filter them by length, if length_range provided
if (!is.null(opt\$length_range)) {

  min_length <- as.integer(strsplit(opt\$length_range, ":")[[1]][1])
  max_length <- as.integer(strsplit(opt\$length_range, ":")[[1]][2])
    
  filtered.ls <- length_filter(data = filtered.ls,
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
psite_offset.dt <- psite(filtered.ls, flanking = 6, extremity = "auto", txt = TRUE, 
                         plot = TRUE, plot_format = "pdf")

# Save offsets for each sample
lapply(unique(psite_offset.dt$sample), export_offsets, df = psite_offset.dt)


# Update reads information according to the inferred P-sites
filtered_psite.ls <- psite_info(filtered.ls, psite_offset.dt, site = "psite",
                                fasta_genome = TRUE, refseq_sep = " ",
                                fastapath = opt\$fasta,
                                gtfpath = gtf)

sample_name.ls <- unique(names(filtered_psite.ls))

# Save psite info for each sample
lapply(sample_name.ls, export_psites, df_list = filtered_psite.ls)

message("Calculating codon and CDS coverage...")

# 1. codon_coverage called by export_codon_coverage_tables computes the number of read footprints or P-sites mapping on each triplet of annotated coding sequences and UTRs. 
# Such data can be exploited to generate occupancy profiles at codon resolution showing the abundance of RPFs along single transcripts
# This function computes transcript-specific codon coverages, defined as the number of either read footprints or P-sites mapping on each triplet of coding sequences and UTRs
lapply(sample_name.ls, export_codon_coverage_tables, df.ls = filtered_psite.ls, annotation.df = annotation.dt)


# 2. Compute the number of P-sites mapping on annotated coding sequences or whole transcripts. 
# Such data can be used as starting point for downstream quantitative analyses (e.g. differential analyses)
# By default, only in-frame P-sites falling in annotated coding sequences are considered 
# no nucleotides at the beginning or at the end of the CDSs are excluded for restricting the analysis to a portion of the original coding sequences. 
# These settings can be modified through the parameters in_frame, start_nts and start_nts. 
# The parameter whole_transcript specifies if whole transcripts should be considered instead of the annotated coding sequence (not used here)

# All over the CDS
cds_coverage_psite.dt <- cds_coverage(filtered_psite.ls, annotation = annotation.dt)

# Compute the number of in-frame P-sites per coding sequences excluding the first 15 codons and the last 10 codons
cds_coverage_psite_window.dt <- cds_coverage(filtered_psite.ls, annotation = annotation.dt, start_nts = opt\$exclude_start, stop_nts = opt\$exclude_stop)
lapply(sample_name.ls, export_cds_coverage_tables, cds_coverage = cds_coverage_psite.dt, cds_window_coverage = cds_coverage_psite_window.dt,
      exclude_start = opt\$exclude_start, exclude_stop = opt\$exclude_stop)


# =========
# Diagnostic plots
# =========

message("Generating diagnostic plots...")

lapply(names(reads.ls), save_length_distribution_plot, dt.ls = reads.ls)

# Metaheatmaps: the abundance of the 5' and 3' extremity of reads mapping on and around the start and the stop codon of annotated CDSs, stratified by their length.
ends_heatmap.gg.ls <- lapply(names(reads.ls), plot_metaheatmap, df_list = reads.ls, annotation = annotation.dt)

# Ribosome profiling data should highlight the CDS of transcripts as the region with the higher percentage of reads. 
lapply(sample_name.ls, save_psite_region_plot, dt.ls = filtered_psite.ls, annotation.df = annotation.dt)

# Plots should show an enrichment of P-sites in the first frame on the coding sequence but not the UTRs, as expected for ribosome protected fragments from protein coding mRNAs.
lapply(sample_name.ls, save_frame_plots, dt.ls = filtered_psite.ls, annotation.df = annotation.dt)

# Trinucleotide periodicity along coding sequences is provided by the function metaprofile_psite. 
# It generates metaprofiles (the merge of single, transcript-specific profiles) based on P-sites mapping around the start and the stop codon of annotated CDSs.
lapply(sample_name.ls, save_metaprofile_psite_plot, df.ls = filtered_psite.ls, annotation.df = annotation.dt)


# Codon usage
lapply(names(filtered_psite.ls), plot_codon_usage, psite_info_ls = filtered_psite.ls)

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
