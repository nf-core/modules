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
    path "versions.yml"                 , emit: versions, topic: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact3: \$(vcontact3 version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/version.json
    touch ${prefix}/stub_database.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcontact3: \$(vcontact3 version)
    END_VERSIONS
    """
}
