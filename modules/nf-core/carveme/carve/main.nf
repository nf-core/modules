process CARVEME_CARVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/carveme:1.6.6--pyhdfd78af_1'
        : 'biocontainers/carveme:1.6.6--pyhdfd78af_1'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.xml")                                                                              , emit: model
    tuple val("${task.process}"), val('carveme'), eval("pip show carveme | sed -n 's/^Version: //p'"), topic: versions, emit: versions_carveme

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    carve \\
        ${fasta} \\
        --output ${prefix}.xml \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.xml
    """
}
