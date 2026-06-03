process GAPSEQ_DOALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/gapseq:2.0.1--5e0dffc1176c5fd2' :
        'quay.io/biocontainers/gapseq:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(medium)

    output:
    tuple val(meta), path("*.RDS")  , emit: model
    tuple val(meta), path("*.tbl")  , emit: tbl
    tuple val(meta), path("*.fna")  , emit: fna      , optional: true
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def medium_arg = medium ? "-m $medium" : ''
    """
    gapseq \\
        doall \\
        -t Bacteria \\
        $medium_arg \\
        -K ${task.cpus} \\
        $args \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_model-filled.RDS
    touch ${prefix}_pathways.tbl
    touch ${prefix}_transporters.tbl
    touch ${prefix}.fna
    touch ${prefix}.log
    """
}
