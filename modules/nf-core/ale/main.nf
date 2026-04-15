process ALE {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ale:20180904--py27ha92aebf_0':
        'biocontainers/ale:20180904--py27ha92aebf_0' }"

    input:
    tuple val(meta), path(asm), path(bam)

    output:
    tuple val(meta), path("*_ALEoutput.txt"), emit: ale
    tuple val("${task.process}"), val('ALE'), val("20180904"), emit: versions_ale, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    ALE \\
        ${args} \\
        ${bam} \\
        ${asm} \\
        ${prefix}_ALEoutput.txt
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ALEoutput.txt
    """

}
