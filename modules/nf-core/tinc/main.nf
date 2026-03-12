#!/usr/bin/env nextflow

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
    template "main_script.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_fit.rds
    touch ${prefix}_plot.rds
    touch ${prefix}.pdf
    touch ${prefix}_qc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtinc: \$(Rscript -e "library(TINC); cat(as.character(packageVersion('TINC')))")
    END_VERSIONS
    """
}
