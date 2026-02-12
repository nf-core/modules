process DSSP_MKDSSP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c9b17fef652145b415471b727cdfd9ed18025bc3551cf296b227f6bafd8a5a1/data':
        'community.wave.seqera.io/library/dssp:4.5.3--019c84afbb8c3190' }"

    input:
    tuple val(meta), path(pdb)
    val(format)

    output:
    tuple val(meta), path("*.{dssp,mmcif}"), emit: dssp
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdssp \\
        $args \\
        --output-format=${format} \\
        ${pdb} \\
        ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dssp: \$(mkdssp --version | sed -n 's/^mkdssp version //p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dssp: \$(mkdssp --version | sed -n 's/^mkdssp version //p')
    END_VERSIONS
    """
}
