process KMCP_COMPUTE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kmcp:0.9.4--h9ee0642_0'
        : 'quay.io/biocontainers/kmcp:0.9.4--h9ee0642_0'}"

    input:
    tuple val(meta), path(sequences, stageAs: "genomes/")

    output:
    tuple val(meta), path("${prefix}"), emit: outdir
    tuple val(meta), path("${prefix}/_info.txt"), emit: info
    tuple val("${task.process}"), val('kmcp'), eval("kmcp version 2>&1 | sed 's/^.*kmcp v//'"), emit: versions_kmcp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmcp \\
        compute \\
        ${args} \\
        --threads ${task.cpus} \\
        --out-dir ${prefix}/ \\
        --in-dir genomes/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/_info.txt
    """
}
