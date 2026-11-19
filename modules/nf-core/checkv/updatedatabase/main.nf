// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CHECKV_UPDATEDATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb71a9ba838eae927761d0c123653534b4a83ef7a9b8fafe3c1de11b5ec71a4/data'
        : 'community.wave.seqera.io/library/checkv_diamond:49f136fbaa08f899'}"

    input:
    tuple val(meta), path (fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/*"), emit: checkv_db
    tuple val("${task.process}"), val("checkv"), eval("checkv -h 2>&1 | sed '1!d;s/^.*CheckV v//;s/:.*//'"), topic: versions, emit: versions_checkv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def checkv_db = db ?: ''
    def update_sequence = fasta ?: ''
    """
    checkv update_database \\
        --threads $task.cpus \\
        $args \\
        $checkv_db \\
        ./$prefix/  \\
        $update_sequence
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/README.txt
    """

}
