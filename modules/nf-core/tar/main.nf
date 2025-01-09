process TAR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/98/98946ea8217c35441352a94f3e0cd1dfa24137c323e8b0f5dfcb3123b465d0b1/data':
        'community.wave.seqera.io/library/bzip2_gzip_lzip_lzop_pruned:5a822ddcf829e7af' }"

    input:
    tuple val(meta), path(input)
    val compress_type

    output:
    tuple val(meta), path("*.tar${compress_type}"), emit: archive
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    valid_compress_types = ['.bz2', '.xz', '.lz', '.lzma', '.lzo', '.zst', '.gz', '']
    if (!compress_type in valid_compress_types) {
        error("ERROR: Invalid compress_type: ${compress_type} for TAR. Set as empty string for no compression. Compression options: ${valid_compress_types.join(", ")}")
    }

    if (compress_type == '.bz2') {
        compress_flag = '--bzip2'
    } else if (compress_type == '.xz') {
        compress_flag = '--xz'
    } else if (compress_type == '.lz') {
        compress_flag = '--lzip'
    } else if (compress_type == '.lzma') {
        compress_flag = '--lzma'
    } else if (compress_type == '.lzo') {
        compress_flag = '--lzop'
    } else if (compress_type == '.zst') {
        compress_flag = '--zstd'
    } else if (compress_type == '.gz') {
        compress_flag = '--gzip'
    } else if (compress_type == '') {
        compress_flag = ''
    } else {
        error("ERROR: Invalid compress_type: ${compress_type} for TAR. Set as empty string for no compression. Compression options: ${valid_compress_types.join(", ")}")
    }

    """
    tar \\
        -c \\
        -h \\
        ${compress_flag} \\
        ${args} \\
        -f ${prefix}.tar${compress_type} \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | grep tar | sed 's/.*) //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip -c > ${prefix}.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | grep tar | sed 's/.*) //g')
    END_VERSIONS
    """
}
