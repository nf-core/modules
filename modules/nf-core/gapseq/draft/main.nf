process GAPSEQ_DRAFT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gapseq:1.4.0--hdfd78af_0' :
        'biocontainers/gapseq:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(pathways), path(transporters)

    output:
    tuple val(meta), path("*.RDS")                                                                 , emit: draft
    tuple val(meta), path("*.log")                                                                 , emit: log    , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq --version | sed "s/gapseq //"')     , emit: versions_gapseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def transporters_arg = transporters ? "-t $transporters" : ''
    """
    gapseq \\
        draft \\
        -r $pathways \\
        $transporters_arg \\
        -c $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-draft.RDS
    """
}
