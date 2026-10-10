process CIRCLATOR_ALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/circlator:1.5.5--py35_0':
       'biocontainers/circlator:1.5.5--py35_0' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(meta2), path(reads)

    output:
    tuple val(meta), path("${prefix}"),        emit: results
    tuple val(meta), path("${prefix}/*.fasta"), emit: fasta
    tuple val("${task.process}"), val('circlator'), eval("circlator version"), topic: versions, emit: versions_circlator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    circlator all \\
        $args \\
        --threads $task.cpus \\
        ${assembly} \\
        ${reads} \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/06.fixstart.fasta
    """
}
