// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_TAXONOMY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/161e2ded0fc1b8d724333808695d7da12a8643668a59b652796e7a51ff9b4f3e/data'
        : 'community.wave.seqera.io/library/mmseqs2:20.ae60e--0746edfc81192f32'}"

    input:
    tuple val(meta), path(db_query)
    path db_target

    output:
    tuple val(meta), path("${prefix}_taxonomy"), emit: db_taxonomy
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    //represents the db_query
    def args3 = task.ext.args3 ?: "*.dbtype"
    //represents the db_target
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_taxonomy

    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_QUERY_PATH_NAME=\$(find -L "${db_query}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_TARGET_PATH_NAME=\$(find -L "${db_target}/" -maxdepth 1 -name "${args3}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        taxonomy \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix}_taxonomy/${prefix} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir -p ${prefix}_taxonomy
    touch ${prefix}_taxonomy/${prefix}.{0..25}
    touch ${prefix}_taxonomy/${prefix}.dbtype
    touch ${prefix}_taxonomy/${prefix}.index

    """
}
