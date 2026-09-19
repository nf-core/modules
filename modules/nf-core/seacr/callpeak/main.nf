process SEACR_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' :
        'quay.io/biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' }"

    input:
    tuple val(meta), path(bedgraph), path(ctrlbedgraph)
    val (threshold)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('scramble'), val('1.3'), topic: versions, emit: versions_scramble
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'")   , topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val('r-base')  , eval("R --version | sed '1!d; s/.*version //; s/ .*//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def function_switch = ctrlbedgraph ? "${ctrlbedgraph}" : "${threshold}"
    """
    SEACR_1.3.sh \\
        ${bedgraph} \\
        ${function_switch} \\
        ${args} \\
        ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
