#!/usr/bin/env Rscript

options(easypar.parallel=FALSE)

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]
print(opt)

# Script #####

library(VIBER)
library(dplyr)
library(tidyr)
library(ggplot2)
library(CNAqc)
library(ggpubr)

samples = substr("$tumour_samples", 2, nchar("$tumour_samples")-1)
samples = strsplit(samples, ", ")[[1]]

cli::cli_text("Patient $meta.patient with samples $tumour_samples.")
cli::cli_text("Input file: $rds_join.")

if ( grepl(".rds\$", tolower("$rds_join")) ) {
    input_obj = readRDS("$rds_join")
    if (class(input_obj) == "m_cnaqc") {
        shared = input_obj %>% get_sample(sample=samples, which_obj="shared")
        joint_table = lapply(names(shared),
                        function(sample_name) {
                          table_s = shared[[sample_name]] %>%
                              CNAqc::subset_by_segment_karyotype("1:1") %>%
                              CNAqc::Mutations() %>%
                              dplyr::mutate(sample_id=sample_name)
                          if (nrow(table_s) == 0) {
                            cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                          }
                          return(table_s)
                        }) %>% dplyr::bind_rows()
        } else {
          cli::cli_alert_warning("Object of class {class(input_obj)} not supported.")
          return()
        }
} else {
  joint_table = read.csv("$rds_join")
}

if (nrow(joint_table) == 0) {
  cli::cli_alert_warning("No samples contain diploid mutations!")
}

cli::cli_text("Joint table created.")

## Read input joint table
input_tab = joint_table %>%
  dplyr::mutate(VAF=replace(VAF, VAF==0, 1e-7))

## Convert the input table into longer format
reads_data = input_tab %>%
  dplyr::select(chr, from, ref, alt, NV, DP, VAF, sample_id,driver_label,is_driver) %>%
  dplyr::rename(gene=driver_label) %>%
  dplyr::rename(driver=is_driver) %>%
  tidyr::pivot_wider(names_from="sample_id",
                     values_from=c("NV","DP","VAF"), names_sep=".",values_fill=0)

## Extract DP (depth)
dp = reads_data %>%
  dplyr::select(dplyr::starts_with("DP")) %>%
  dplyr::rename_all(function(x) stringr::str_remove_all(x,"DP."))

## Extract NV (alt_counts)
nv = reads_data %>%
  dplyr::select(dplyr::starts_with("NV")) %>%
  dplyr::rename_all(function(x) stringr::str_remove_all(x,"NV."))

# Standard fit
viber_K = as.integer(opt[["K"]])

cli::cli_text("Starting standard fit")
st_fit = VIBER::variational_fit(nv, dp,
                                K=viber_K,
                                data=reads_data)
st_fit[["description"]]="$meta.patient"
cli::cli_text("End standard fit")
best_fit = best_fit_heuristic = st_fit

# If all clusters are removed -> keep the origianl best fit
tryCatch(expr = {
    # Apply the heuristic
    best_fit_heuristic = VIBER::choose_clusters(st_fit,
        binomial_cutoff=as.numeric(opt[["binomial_cutoff"]]),
        dimensions_cutoff=as.integer(opt[["dimensions_cutoff"]]),
        pi_cutoff=as.numeric(opt[["pi_cutoff"]]),
        re_assign=as.logical(opt[["re_assign"]])
        )
    }, error = function(e) {
        print(e)
        best_fit_heuristic <<- st_fit
    } )

# Save fits
saveRDS(best_fit, file=paste0(opt[["prefix"]], "_viber_best_fit.rds"))
saveRDS(best_fit_heuristic, file = paste0(opt[["prefix"]], "_viber_best_heuristic_fit.rds"))

# Save plots
n_samples = ncol(best_fit[["x"]]) - 1
if (n_samples >1) {
  print("multisample mode on")
  plot_fit = plot(best_fit)
  plot_fit_heuristic = plot(best_fit_heuristic)
} else {
  plot_fit = plot_mixing_proportions(best_fit)
  plot_fit_heuristic = plot_mixing_proportions(best_fit_heuristic)
}

saveRDS(plot_fit, file=paste0(opt[["prefix"]], "_viber_best_fit_plots.rds"))
saveRDS(plot_fit_heuristic, file=paste0(opt[["prefix"]], "_viber_best_heuristic_fit_plots.rds"))

# Save report plot
n_samples = ncol(best_fit[["x"]]) - 1
marginals = multivariate = ggplot()

try(expr = {marginals <<- VIBER::plot_1D(best_fit)} )
try(expr = {multivariate = plot(best_fit)})
try(expr = {multivariate = ggpubr::ggarrange(plotlist = multivariate)})
top_p = ggpubr::ggarrange(plotlist = list(marginals, multivariate), widths=ifelse(n_samples>2, c(1,2), c(2,1)))

mix_p = VIBER::plot_mixing_proportions(best_fit)
binom = VIBER::plot_peaks(best_fit)
elbo = VIBER::plot_ELBO(best_fit)
bottom_p = ggpubr::ggarrange(plotlist = list(mix_p, binom, elbo), widths = c(1,2,1), nrow = 1)

report_fig = ggpubr::ggarrange(top_p, bottom_p, nrow = 2, heights = ifelse(n_samples>2, c(3,1), c(2,1)))
saveRDS(report_fig, file=paste0(opt[["prefix"]], "_viber_report.rds"))
ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_viber_report.pdf"), height=210, width=210, units="mm", dpi = 200)
ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_viber_report.png"), height=210, width=210, units="mm", dpi = 200)

# version export
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
viber_version = sessionInfo()\$otherPkgs\$VIBER\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
tidyr_version = sessionInfo()\$otherPkgs\$tidyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
ggpubr_version = sessionInfo()\$otherPkgs\$ggpubr\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    VIBER:", viber_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    tidyr:", tidyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    ggpubr:", ggpubr_version), f)
close(f)
