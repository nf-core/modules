process TINC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tinc%3A0.1.0--r44hdfd78af_0':
        'biocontainers/r-tinc:0.1.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(cna_rds), path(vcf_rds)

    output:
    tuple val(meta), path("*_fit.rds"),     emit: rds
    tuple val(meta), path("*_plot.rds"),    emit: plot_rds
    tuple val(meta), path("*.pdf"),         emit: plot_pdf
    tuple val(meta), path("*_qc.csv"),      emit: tinc_csv
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                                = task.ext.args     ?: ''                                   
    def prefix                              = task.ext.prefix   ?: "${meta.id}"
    def normal_contamination_level          = args!='' && args.normal_contamination_level           ?  "$args.normal_contamination_level" : ""

    """
    #!/usr/bin/env Rscript

    require(CNAqc)
    require(tidyverse)
    require(TINC)

    all_mutations = readRDS("$vcf_rds")
    samples = names(all_mutations)

    tumor_sample = "$meta.tumour_sample"
    normal_sample = "$meta.normal_sample"

    tumor_mutations = all_mutations[[tumor_sample]]\$mutations %>%
        select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
        filter(!is.na(DP)) %>%
        rename(t_alt_count = NV, t_ref_count = NR, t_tot_count = DP, t_vaf = VAF)

    normal_mutations = all_mutations[[normal_sample]]\$mutations %>%
        select(chr, from, to, ref, alt, NV, DP, NR, VAF) %>%
        filter(!is.na(DP)) %>%
        rename(n_alt_count = NV, n_ref_count = NR, n_tot_count = DP, n_vaf = VAF)

    input_mut = dplyr::full_join(tumor_mutations, normal_mutations, by = c("chr", "from", "to", "ref", "alt")) %>%
        mutate(t_vaf = case_when(is.na(t_vaf) ~ 1e-5, .default = t_vaf)) %>%
        mutate(n_vaf = case_when(is.na(n_vaf) ~ 1e-5, .default = n_vaf)) %>%
        mutate(t_vaf = as.numeric(t_vaf), n_vaf = as.numeric(n_vaf)) %>%
        filter(!(is.na(t_alt_count))) %>%
        filter(!(is.na(n_alt_count))) %>%
        #mutate(t_alt_count = case_when(is.na(t_alt_count) ~ 0, .default = t_alt_count)) %>%
        #mutate(t_ref_count = case_when(is.na(t_ref_count) ~ 0, .default = t_ref_count)) %>%
        #mutate(n_alt_count = case_when(is.na(n_alt_count) ~ 0, .default = n_alt_count)) %>%
        #mutate(n_ref_count = case_when(is.na(n_ref_count) ~ 0, .default = n_ref_count)) %>%
        mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count), n_tot_count = as.numeric(n_tot_count), n_ref_count = as.numeric(n_ref_count)) %>%
        filter(t_vaf > 0)

    CNAs = readRDS("$cna_rds")\$segments
    TINC_fit = TINC::autofit(input = input_mut,
                    cna = CNAs)

    tinc_plot = plot(TINC_fit)

    qc_res = TINC:::classification_normal(TINC_fit\$TIN)

    if (qc_res[["level"]] >= eval(parse(text="$normal_contamination_level"))) {
        sample_contamination = tibble(
            sample = tumor_sample,
            normal_contamination = qc_res[["level"]],
            normal_contamination_flag = 1
        )
    } else {
        sample_contamination = tibble(sample = tumor_sample,
        normal_contamination = qc_res[["level"]],
        normal_contamination_flag = 0
        )
    }

    write.table(file = paste0("${prefix}", "_qc.csv"), sep = ",", x = sample_contamination, col.names = T, row.names = F, quote = F)

    saveRDS(file = paste0("${prefix}", "_plot.rds"), object = tinc_plot)
    ggplot2::ggsave(plot = tinc_plot, filename = paste0("${prefix}", "_plot.pdf"), width = 210, height = 297, units="mm", dpi = 200)
    saveRDS(file = paste0("${prefix}", "_fit.rds"), object = TINC_fit)

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

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    
    touch ${prefix}_fit.rds
    touch ${prefix}_plot.rds
    touch ${prefix}.pdf
    touch ${prefix}_qc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-r-tinc:
    END_VERSIONS
    """
}
