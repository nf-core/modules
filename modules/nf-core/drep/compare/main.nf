process DREP_COMPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/drep:3.6.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/drep:3.6.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${prefix}"), emit: directory
    tuple val("${task.process}"), val("drep"), eval("dRep | head -n 2 | sed 's/.*v//g;s/ .*//g' | tail -n 1"), emit:versions_drep, topic:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    dRep \\
        compare \\
        ${prefix} \\
        -p ${task.cpus} \\
        ${args} \\
        -g fastas/*
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p ${prefix}
    """
}
