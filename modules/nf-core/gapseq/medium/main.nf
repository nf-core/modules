process GAPSEQ_MEDIUM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/933e301b11c1ec1699da6382e9e35b0e4e31edb80763eb2fa1b69ad7d6d1e5c7/data'
:         'community.wave.seqera.io/library/gapseq:2.1.0--c32b876ebb5e5f5b' }"

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
        --model $draft \
        --pathway.pred $pathways \
        $args \
        --output.file ${prefix}-medium.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-medium.csv
    """
}
