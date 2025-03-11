process CHECKM2_DATABASEDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aria2:1.36.0':
        'biocontainers/aria2:1.36.0' }"

    input:
    val(db_zenodo_id)

    output:
    tuple val(meta), path("checkm2_db_v${db_version}.dmnd"), emit: database
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    zenodo_id       = db_zenodo_id ?: 5571251  // Default to latest version if no ID provided
    api_data        = (new groovy.json.JsonSlurper()).parseText(file("https://zenodo.org/api/records/${zenodo_id}").text)
    db_version      = api_data.metadata.version
    checksum        = api_data.files[0].checksum.replaceFirst(/^md5:/, "md5=")
    meta            = [id: 'checkm2_db', version: db_version]
    """
    # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
    # So it's necessary to download the database manually
    aria2c \
        ${args} \
        --checksum ${checksum} \
        https://zenodo.org/records/${zenodo_id}/files/checkm2_database.tar.gz

    tar -xzf checkm2_database.tar.gz
    db_path=\$(find -name *.dmnd)
    mv \$db_path checkm2_db_v${db_version}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """

    stub:
    """
    touch checkm_db.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
