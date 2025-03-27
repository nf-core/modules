process CHECKM2_DATABASEDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data':
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
    val(db_zenodo_id)

    output:
    tuple val(meta), path("checkm2_db_v${db_version}.dmnd"), emit: database
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    zenodo_id       = db_zenodo_id ?: 14897628  // Default to version 3 if no ID provided
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
    db_version = 0
    meta       = [id: 'checkm2_db', version: db_version]
    """
    touch checkm2_db_v${db_version}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
