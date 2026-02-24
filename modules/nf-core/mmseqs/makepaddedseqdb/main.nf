process MMSEQS_MAKEPADDEDSEQDB {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    tuple val(meta), path(prefix)

    output:
    tuple val(meta), path("${padded_prefix}/"), emit: db_padded
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    padded_prefix = task.ext.prefix ?: "${meta.id}_gpu"
    """
    mkdir -p ${padded_prefix}
    mmseqs \\
        makepaddedseqdb \\
        ${prefix}/${prefix} \\
        ${padded_prefix}/${padded_prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    padded_prefix = task.ext.prefix ?: "${meta.id}_gpu"
    """
    echo ${args}
    mkdir -p ${padded_prefix}
    touch ${padded_prefix}/${padded_prefix}
    touch ${padded_prefix}/${padded_prefix}.dbtype
    touch ${padded_prefix}/${padded_prefix}.index
    touch ${padded_prefix}/${padded_prefix}.lookup
    touch ${padded_prefix}/${padded_prefix}_h
    touch ${padded_prefix}/${padded_prefix}_h.dbtype
    touch ${padded_prefix}/${padded_prefix}_h.index
    """
}
