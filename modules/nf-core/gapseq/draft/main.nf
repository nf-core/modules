process GAPSEQ_DRAFT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/gapseq:2.0.1--5e0dffc1176c5fd2' :
        'quay.io/biocontainers/gapseq:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reactions), path(transporters), path(pathways)

    output:
    tuple val(meta), path("*.RDS")  , emit: draft
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def transporters_arg = transporters ? "-t $transporters" : ''
    """
    gapseq \\
        draft \\
        -r $reactions \\
        $transporters_arg \\
        -p $pathways \\
        -n $prefix \
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-draft.RDS
    """
}
