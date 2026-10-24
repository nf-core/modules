// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GCTA_ADDGRMS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56e7d9dd85fa648387dcb28e6f9c2eb293c5b2f444f5419266be7732cac23f3d/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--e22d27187f066d3c'}"

    input:
    tuple val(meta), path(mgrm_file), path(grm_files)

    output:
    tuple val(meta), path("*.grm.*"), emit: combined_grm
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gcta \\
        --mgrm ${mgrm_file} \\
        --make-grm \\
        --out ${prefix} \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grm.id
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.N.bin
    """
}
