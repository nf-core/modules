process HMMER_HMMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b171806a7b780537073c1e9ffe50885ba7917d8072b5941a125698df1253c4f7/data'
:         'community.wave.seqera.io/library/hmmer:3.4--ab3efa7c16f46920' }"

    input:
    tuple val(meta), path(hmmfile)

    output:
    tuple val(meta), path("*.h3?"), emit: compressed_db
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | sed '2!d;s/^# HMMER *//;s/ .*//'"), emit: versions_hmmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    hmmpress \\
        $args \\
        ${hmmfile}
    """

    stub:
    def prefix = task.ext.prefix ?: "stub"

    """
    touch ${prefix}.h3m
    touch ${prefix}.h3i
    touch ${prefix}.h3f
    touch ${prefix}.h3p
    """
}
