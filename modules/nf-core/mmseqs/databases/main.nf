// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_DATABASES {
    tag "${database}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7df92cbef82c41ca1fcdc2959e9beeeeb6106099accfdaaed468b2990587061/data'
        : 'community.wave.seqera.io/library/mmseqs2_wget:4ce0f9394a204754'}"

    input:
    val database

    output:
    path "${prefix}/", emit: database
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir ${prefix}/

    mmseqs databases \\
        ${database} \\
        ${prefix}/database \\
        tmp/ \\
        --threads ${task.cpus} \\
        ${args}

    """

    stub:
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir ${prefix}/

    touch ${prefix}/database
    touch ${prefix}/database.dbtype
    touch ${prefix}/database_h
    touch ${prefix}/database_h.dbtype
    touch ${prefix}/database_h.index
    touch ${prefix}/database.index
    touch ${prefix}/database.lookup
    touch ${prefix}/database_mapping
    touch ${prefix}/database.source
    touch ${prefix}/database_taxonomy
    touch ${prefix}/database.version

    """
}
