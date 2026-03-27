process DSSP_MKDSSP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e45865baee2b1563b471eb47b62f04ad6e49a862bcb815f8944b5bdba8cc33e1/data':
        'community.wave.seqera.io/library/dssp:4.5.8--9dc7af262c1d3dd7' }"

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
