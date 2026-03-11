process CNAQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-cnaqc%3A1.1.3--r44hdfd78af_0':
        'biocontainers/r-cnaqc:1.1.3--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(snv_rds), path(cna_rds), val(tumour_sample)

    output:
    tuple val(meta), path("*_qc.rds"),                                  emit: qc_rds
    tuple val(meta), path("*_data_plot.rds"), path("*_qc_plot.rds"),    emit: data_plot_rds
    tuple val(meta), path("*_qc_plot.rds"),                             emit: qc_plot_rds
    tuple val(meta), path("*_data.pdf"),                                emit: plot_pdf_data
    tuple val(meta), path("*_qc.pdf"),                                  emit: plot_pdf_qc
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    "echo ${args}"

    template "main_script.R"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ''
    """
    touch ${prefix}_qc.rds
    touch ${prefix}_data_plot.rds
    touch ${prefix}_qc_plot.rds
    touch ${prefix}_qc.pdf
    touch ${prefix}_data.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnaqc: \$(Rscript -e "library(CNAqc); cat(as.character(packageVersion('CNAqc')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
