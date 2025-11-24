#!/usr/bin/env Rscript

# parse_args = function(x) {
#     x = gsub("\\\\[","",x)
#     x = gsub("\\\\]","",x)
#     # giving errors when we have lists like c(xxx, xxx) since it will separate it
#     # args_list = unlist(strsplit(x, ', ')[[1]])
#     args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
#     # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
#     args_vals = lapply(args_list, function(x) {
#         x_splt = strsplit(x, split=":")[[1]]
#         c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
#     })

#     # Ensure the option vectors are length 2 (key/ value) to catch empty ones
#     args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

#     parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
#     parsed_args[! is.na(parsed_args)]
# }

# opt = list(
#     prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
# )
# args_opt = parse_args('$task.ext.args')
# for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]
# print(opt)

# parse arguments
parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

# Set defaults and classes

opt <- list(
  prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),
  genome = 'GRCh38',
  genome_coords = NULL,
  karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2"),
  min_absolute_karyotype_mutations = 100,
  matching_strategy = 'rightmost',
  min_karyotype_size = 0,
  p_binsize_peaks = 0.005,
  matching_epsilon = NULL,
  purity_error = 0.05,
  VAF_tolerance = 0.015,
  n_bootstrap = 1,
  kernel_adjust = 1,
  KDE = TRUE,
  starting_state_subclonal_evolution = "1:1",
  cluster_subclonal_CCF = FALSE,
  min_VAF = 0,
  muts_per_karyotype = 25,
  cutoff_QC_PASS = 0.1,
  method = "ENTROPY"
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{

    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}

# Script #####

# load libraries

library(dplyr)
library(CNAqc)

SNV = readRDS("$snv_rds") %>%
  purrr::pluck("$tumour_sample", "mutations") %>%
  dplyr::mutate(mutation_id = paste(chr,from,to,ref,alt,sep = ':'))

CNA = readRDS("$cna_rds")

x = CNAqc::init(
  mutations = SNV,
  cna = CNA\$segments,
  purity = CNA\$purity ,
  ref = opt[["genome"]])

x = CNAqc::analyze_peaks(x,
                          matching_strategy = opt[["matching_strategy"]],
                          min_absolute_karyotype_mutations = as.numeric(opt[["min_absolute_karyotype_mutations"]]),
                          purity_error = as.numeric(opt[["purity_error"]]))

x = CNAqc::compute_CCF(x,
                        muts_per_karyotype = as.numeric(opt[["muts_per_karyotype"]])
)

# this is needed in order to plot the results without the 0 VAF mutations
tmp_x <- x
mut <- CNAqc::Mutations(tmp_x) %>%
  dplyr::filter(VAF > 0)
tmp_x\$mutations <- mut

pl = ggpubr::ggarrange(
  CNAqc::plot_data_histogram(tmp_x, which = 'VAF'),
  CNAqc::plot_data_histogram(tmp_x, which = 'DP'),
  CNAqc::plot_data_histogram(tmp_x, which = 'NV'),
  CNAqc::plot_data_histogram(tmp_x, which = 'CCF'),
  ncol = 2,
  nrow = 2
)

pl_exp = ggpubr::ggarrange(
  plotlist = list(CNAqc::plot_gw_counts(tmp_x),
                  CNAqc::plot_gw_vaf(tmp_x, N = 10000),
                  CNAqc::plot_gw_depth(tmp_x, N = 10000),
                  CNAqc::plot_segments(tmp_x),
                  pl),
  nrow = 5,
  heights = c(.5,.5,.5,1,5)
)
pl_exp = ggpubr::annotate_figure(pl_exp, top = ggpubr::text_grob("$tumour_sample", size = 14))

pl_qc = ggpubr::ggarrange(
  plotlist = list(
    CNAqc::plot_peaks_analysis(tmp_x, what = 'common', empty_plot = FALSE),
    CNAqc::plot_qc(tmp_x),
    CNAqc::plot_CCF(tmp_x, assembly_plot = TRUE, empty_plot = FALSE)),
  nrow = 3,
  heights = c(1,1.5,1))
pl_qc = ggpubr::annotate_figure(pl_qc, top = ggpubr::text_grob("$tumour_sample", size = 14))

saveRDS(object = x, file = paste0(opt[["prefix"]], "_qc.rds"))
saveRDS(object = pl_exp, file = paste0(opt[["prefix"]], "_data_plot.rds"))
saveRDS(object = pl_qc, file = paste0(opt[["prefix"]], "_qc_plot.rds"))

ggplot2::ggsave(plot = pl_exp, filename = paste0(opt[["prefix"]], "_data.pdf"), width = 210, height = 297, units="mm", dpi = 200)
ggplot2::ggsave(plot = pl_qc, filename = paste0(opt[["prefix"]], "_qc.pdf"), width = 210, height = 297, units="mm", dpi = 200)

# version export
f <- file("versions.yml","w")
dplyr_version <- sessionInfo()\$otherPkgs\$dplyr\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
close(f)
