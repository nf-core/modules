process VCONTACT3_PREPAREDDATABASES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcontact3:3.1.6--pyhdfd78af_1':
        'quay.io/biocontainers/vcontact3:3.1.6--pyhdfd78af_1' }"

    input:
    tuple val(meta), val(db_version)

    output:
    tuple val(meta), path ("${prefix}/"), emit: database
    tuple val("${task.process}"), val('vcontact3'), eval('vcontact3 version'), emit: versions_vcontact3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    vcontact3 prepare_databases \\
        --get-version ${db_version} \\
        --set-location ${prefix}  \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/version.json
    touch ${prefix}/stub_database.db
    """
}
