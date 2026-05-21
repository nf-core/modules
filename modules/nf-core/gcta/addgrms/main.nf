process GCTA_ADDGRMS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(mgrm_file), path(grm_files)

    output:
    tuple val(meta), path("*.grm.*"), emit: combined_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    gcta \\
        --mgrm ${mgrm_file} \\
        --make-grm \\
        --out ${prefix} \\
        --thread-num ${task.cpus} \\
        ${extra_args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grm.id
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.N.bin
    """
}
