process GCTA_KEEP {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(grm_files)
    tuple val(meta2), path(keep_file)

    output:
    tuple val(meta), path("${prefix}.grm.*"), emit: filtered_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_keep"
    """
    gcta \\
        --grm ${meta.id} \\
        --keep ${keep_file} \\
        --make-grm \\
        --out ${prefix} \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_keep"
    """
    touch ${prefix}.grm.id
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.N.bin
    """
}
