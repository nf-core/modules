process GATK4_ANALYZECOVARIATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/3498e3d67f977e18255b0b71560dd548885309f62b736dd2c782b131cf931008/data'
        : 'community.wave.seqera.io/library/gatk4-lite_r-base_r-ggplot2_r-gplots_r-gsalib:aad8b6b5e895ebef'}"

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
