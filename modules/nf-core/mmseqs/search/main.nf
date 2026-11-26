// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_SEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/161e2ded0fc1b8d724333808695d7da12a8643668a59b652796e7a51ff9b4f3e/data'
        : 'community.wave.seqera.io/library/mmseqs2:20.ae60e--0746edfc81192f32'}"

    input:
    tuple val(meta), path(db_query)
    tuple val(meta2), path(db_target)

    output:
    tuple val(meta), path("${prefix}/"), emit: db_search
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    def args3 = task.ext.args3 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"

    if ("${db_query}" == "${prefix}" || "${db_target}" == "${prefix}") {
        error("Input and output names of databases are the same, set prefix in module configuration to disambiguate!")
    }
    """
    mkdir -p ${prefix}

    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_QUERY_PATH_NAME=\$(find -L "${db_query}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_TARGET_PATH_NAME=\$(find -L "${db_target}/" -maxdepth 1 -name "${args3}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        search \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix}/${prefix} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${db_query}" == "${prefix}" || "${db_target}" == "${prefix}") {
        error("Input and output names of databases are the same, set prefix in module configuration to disambiguate!")
    }
    """
    echo ${args}
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.{0..9}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index

    """
}
