def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This module is no longer recommended for use. It is recommended to use nf-core/modules/untar
"""
process UNTARFILES {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/88e03525287eaeb8bb74114aaee2c67118c1cdcfb99ee52e3ddc71a1acce35d4/data' :
        'community.wave.seqera.io/library/grep_sed_tar:db2951cd23a1ffde' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}/**") , emit: files
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert true: deprecation_message
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName.toString().replaceFirst(/\.tar$/, ""))

    """
    mkdir $prefix

    tar \\
        -C $prefix \\
        -xavf \\
        $args \\
        $archive \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    assert true: deprecation_message
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch ${prefix}/file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
