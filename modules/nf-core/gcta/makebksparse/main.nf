process GCTA_MAKEBKSPARSE {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(grm_files)
    val cutoff

    output:
    tuple val(meta), path("*_sp.grm.*"), emit: sparse_grm_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dense_prefix = meta.id
    def sparse_prefix = "${prefix}_sp"
    """
    gcta \\
        --grm ${dense_prefix} \\
        --make-bK-sparse ${cutoff} \\
        --out ${sparse_prefix} \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sparse_prefix = "${prefix}_sp"
    """
    touch ${sparse_prefix}.grm.id
    touch ${sparse_prefix}.grm.sp
    """
}
