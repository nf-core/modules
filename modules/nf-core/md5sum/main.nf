process MD5SUM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(files)
    val as_separate_files

    output:
    tuple val(meta), path("*.md5"), emit: checksum
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // will only use when as_separate_files = false
    if (as_separate_files) {
        """
        find -L * -maxdepth 0 -type f \\
            ! -name '*.md5' \\
            -exec sh -c 'md5sum ${args} "\$1" > "\$1.md5"' _ "{}" \\;

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            md5sum: \$( md5sum --version | sed '1!d; s/.* //' )
        END_VERSIONS
        """
    }
    else {
        """
        find -L * -type f \\
            ! -name '*.md5' \\
            -exec md5sum ${args} "{}" + \\
            > ${prefix}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            md5sum: \$( md5sum --version | sed '1!d; s/.* //' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (as_separate_files) {
        """
        find -L * -type f \\
            ! -name '*.md5' \\
            -exec sh -c 'touch "\$1.md5"' _ "{}" \\;

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            md5sum: \$( md5sum --version | sed '1!d; s/.* //' )
        END_VERSIONS
        """
    }
    else {
        """
        touch ${prefix}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            md5sum: \$( md5sum --version | sed '1!d; s/.* //' )
        END_VERSIONS
        """
    }
}
