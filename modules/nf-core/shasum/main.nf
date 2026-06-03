process SHASUM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.sha256"), emit: checksum
    tuple val("${task.process}"), val('sha256sum'), eval("sha256sum --version 2>&1 | head -n 1 | sed 's/.* //'"), emit: versions_sha256sum, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sha256sum \\
        ${args} \\
        ${file} \\
        > ${file}.sha256
    """

    stub:
    """
    touch ${file}.sha256
    """
}
