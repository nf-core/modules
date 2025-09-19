#!/usr/bin/env Rscript

require(CNAqc)
require(tidyverse)
require(dplyr)
require(TINC)

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
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

all_mutations = readRDS("$vcf_rds")
samples = names(all_mutations)

tumor_sample = "$meta.tumour_sample"
normal_sample = "$meta.normal_sample"

tumor_mutations = all_mutations[[tumor_sample]][['mutations']] %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
    dplyr::filter(!is.na(DP)) %>%
    dplyr::rename(t_alt_count = NV, t_ref_count = NR, t_tot_count = DP, t_vaf = VAF)

normal_mutations = all_mutations[[normal_sample]][['mutations']] %>%
    dplyr::select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
    dplyr::filter(!is.na(DP)) %>%
    dplyr::rename(n_alt_count = NV, n_ref_count = NR, n_tot_count = DP, n_vaf = VAF)

input_mut = dplyr::full_join(tumor_mutations, normal_mutations, by = c("chr", "from", "to", "ref", "alt")) %>%
    dplyr::mutate(t_vaf = case_when(is.na(t_vaf) ~ 1e-5, .default = t_vaf)) %>%
    dplyr::mutate(n_vaf = case_when(is.na(n_vaf) ~ 1e-5, .default = n_vaf)) %>%
    dplyr::mutate(t_vaf = as.numeric(t_vaf), n_vaf = as.numeric(n_vaf)) %>%
    dplyr::filter(!(is.na(t_alt_count))) %>%
    dplyr::filter(!(is.na(n_alt_count))) %>%
    dplyr::mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count), n_tot_count = as.numeric(n_tot_count), n_ref_count = as.numeric(n_ref_count)) %>%
    dplyr::filter((n_ref_count + n_alt_count) > 0) %>%
    dplyr::filter(t_vaf > 0)

CNAs = readRDS("$cna_rds")[['segments']]
TINC_fit = TINC::autofit(input = input_mut,
                cna = CNAs)

tinc_plot = plot(TINC_fit)

qc_res = TINC:::classification_normal(TINC_fit[['TIN']])

if (qc_res[["level"]] >= eval(parse(text=opt[["normal_contamination_level"]]))) {
    sample_contamination = dplyr::tibble(
        sample = tumor_sample,
        normal_contamination = qc_res[["level"]],
        normal_contamination_flag = 1
    )
} else {
    sample_contamination = dplyr::tibble(sample = tumor_sample,
    normal_contamination = qc_res[["level"]],
    normal_contamination_flag = 0
    )
}

write.table(file = paste0(opt[["prefix"]], "_qc.csv"), sep = ",", x = sample_contamination, col.names = T, row.names = F, quote = F)

saveRDS(file = paste0(opt[["prefix"]], "_plot.rds"), object = tinc_plot)
ggplot2::ggsave(plot = tinc_plot, filename = paste0(opt[["prefix"]], "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
saveRDS(file = paste0(opt[["prefix"]], "_fit.rds"), object = TINC_fit)

# version export
f <- file("versions.yml","w")
tidyverse_version <- sessionInfo()\$otherPkgs\$tidyverse\$Version
cnaqc_version <- sessionInfo()\$otherPkgs\$CNAqc\$Version
tinc_version <- sessionInfo()\$otherPkgs\$TINC\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    tidyverse:", tidyverse_version), f)
writeLines(paste("    TINC:", tinc_version), f)
close(f)
