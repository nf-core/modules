process HMMER_HMMSTAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b171806a7b780537073c1e9ffe50885ba7917d8072b5941a125698df1253c4f7/data'
:         'community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920' }"

    input:
    tuple val(meta), path(hmm)

    output:
    tuple val(meta), path("${prefix}.hmmstat.txt"), emit: stats
    tuple val("${task.process}"), val('hmmer'), eval("hmmstat -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), topic: versions, emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hmmstat \\
        $args \\
        $hmm \\
        > ${prefix}.hmmstat.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hmmstat.txt
    """
}
