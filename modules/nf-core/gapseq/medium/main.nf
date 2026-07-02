process GAPSEQ_MEDIUM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/gapseq:2.1.0--31c8824b3592beaf' :
        'quay.io/biocontainers/gapseq:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(draft), path(pathways)

    output:
    tuple val(meta), path("*.csv") , emit: medium
    tuple val(meta), path("*.log") , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gapseq \
        medium \
        -m $draft \
        -p $pathways \
        $args \
        -o ${prefix}-medium.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-medium.csv
    """
}
