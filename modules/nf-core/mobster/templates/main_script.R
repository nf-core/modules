#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
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


# Script #####

library(CNAqc)
library(mobster)
library(dplyr)
library(ggplot2)

description = "$meta.patient"
samples = strsplit(x="$meta.tumour_sample", ",") %>% unlist()  # list of samples

## read mCNAqc object
if ( grepl(".rds\$", tolower("$rds_join")) ) {
    obj = readRDS("$rds_join")
    if (class(obj) == "m_cnaqc") {
        original = obj %>% get_sample(sample=samples, which_obj="original")
        input_table = lapply(names(original),
                             function(sample_name) {
                                 purity = original[[sample_name]][["purity"]]
                                 original[[sample_name]] %>%
                                     # keep only mutations on the diploid karyotype
                                     CNAqc::subset_by_segment_karyotype("1:1") %>%
                                     CNAqc::Mutations() %>%
                                     dplyr::mutate(sample_id=sample_name, purity=purity)
                                 }) %>% dplyr::bind_rows()
    } else {
        cli::cli_abort("Object of class {class($rds_join)} not supported.")
    }
} else {
    input_table = read.csv("$rds_join")
}

## Function to run a single MOBSTER fit
run_mobster_fit = function(joint_table, descr) {
    purity = unique(joint_table[["purity"]]) %>% as.numeric()

    if(length(purity) > 1)
        cli::cli_abort("Multiple values of purity are found in one sample (sample ID: {descr}).")

    # remove mutations with VAF < min_VAF
    inp_tb = joint_table %>%
        dplyr::filter(VAF < 1) %>%
        dplyr::mutate(adj_VAF=VAF/!!purity) %>%
        dplyr::filter(adj_VAF>=as.numeric(opt[["min_VAF"]])) %>%
        dplyr::filter(karyotype=="1:1") %>%
        dplyr::select(-adj_VAF)

    mobster_fit(x = inp_tb,
                K = eval(parse(text=opt[["K"]])),
                tail = eval(parse(text=opt[["tail"]])),
                pi_cutoff = as.numeric(opt[["pi_cutoff"]]),
                N_cutoff = as.integer(opt[["n_cutoff"]]),
                parallel = FALSE,
                description = descr)
}

## Run MOBSTER on each sample
lapply(samples, function(sample_name) {

    fit = run_mobster_fit(joint_table=input_table %>% dplyr::filter(sample_id == !!sample_name),
                          descr=paste(description, sample_name, sep=":"))

    if (any(fit[["fits.table"]][["tail"]])) {
        evoparams = evolutionary_parameters(fit)
        fit[["evolutionary_parameters"]] = evoparams
        fit[["best"]][["evolutionary_parameters"]] = evoparams
    }

    best_fit = fit[["best"]]
    plot_fit = plot(best_fit)

    saveRDS(object=fit, file=paste0(opt[["prefix"]], "_mobster_fit.rds"))
    saveRDS(object=best_fit, file=paste0(opt[["prefix"]], "_mobster_best_fit.rds"))
    saveRDS(object=plot_fit, file=paste0(opt[["prefix"]], "_mobster_best_plots.rds"))

    # save report plots
    report_fig = mobster::plot_model_selection(fit)
    saveRDS(report_fig, file=paste0(opt[["prefix"]], "_mobster_report.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.pdf"), height=210, width=210, units="mm", dpi = 200)
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.png"), height=210, width=210, units="mm", dpi = 200)
})


# Version export #####
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    mobster:", mobster_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
close(f)
