// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MMSEQS_TSV2EXPROFILEDB {
    tag "${database}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/161e2ded0fc1b8d724333808695d7da12a8643668a59b652796e7a51ff9b4f3e/data'
        : 'community.wave.seqera.io/library/mmseqs2:20.ae60e--0746edfc81192f32'}"

    input:
    path database

    output:
    path (database), emit: db_exprofile
    tuple val("${task.process}"), val('mmseqs'), eval('mmseqs version'), topic: versions, emit: versions_mmseqs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    DB_PATH_NAME=\$(find -L "${database}/" -name "*_seq.tsv" | sed 's/_seq\\.tsv\$//')

    mmseqs tsv2exprofiledb \\
        \${DB_PATH_NAME} \\
        "\${DB_PATH_NAME}_db" \\
        ${args}

    """

    stub:
    """
    DB_PATH_NAME=\$(find -L "${database}/" -name "*_seq.tsv" | sed 's/_seq\\.tsv\$//')
    touch "\${DB_PATH_NAME}_db"

    """
}
