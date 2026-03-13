process COWPY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1a1f6ff027a12fc86ad7bb9f8781c06b6fd7c8b81a9ecde90cda09335927b0fe/data'
        : 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'}"

    input:
    tuple val(meta), path(text)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt
    tuple val("${task.process}"), val('cowpy'), val("1.1.5"), emit: versions_cowpy, topic: versions
    // WARN: Version information not provided by tool on CLI. Plaease update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${text} | cowpy ${args} > ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
