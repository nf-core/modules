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
        'oras://community.wave.seqera.io/library/grep_sed_tar:40b34489f8e82876' :
        'community.wave.seqera.io/library/grep_sed_tar:16f6591cd62505b3' }"

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
