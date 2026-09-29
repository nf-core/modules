process GATK4_ANALYZECOVARIATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/57141fe65436d371ade52d3700cc6718ab80b6ddd8a1bc8c6269c3e84b387a41/data'
        : 'community.wave.seqera.io/library/gatk4-main_r-base_r-ggplot2_r-gplots_r-gsalib:dd924f4a430b36fc'}"

    input:
    tuple val(meta), path(before_table), path(after_table), path(additional_table)

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.csv"), emit: data
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def third_table = additional_table ? "-bqsr ${additional_table}" : ""

    """
    gatk AnalyzeCovariates \\
        -before ${before_table} \\
        -after ${after_table} \\
        ${third_table} \\
        -csv ${meta.id}.csv \\
        -plots ${meta.id}.pdf \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''

    """
    echo ${args}

    touch ${meta.id}.csv
    touch ${meta.id}.pdf
    """
}
