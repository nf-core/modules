import groovy.json.JsonSlurper

process CHECKM2_DATABASEDOWNLOAD {
    label 'process_single'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    output:
    tuple val(meta), path("checkm2_db_v${db_version}.dmnd"), emit: database
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    zenodo_id = 5571251
    def jsonSlurper = new JsonSlurper()
    db_version = jsonSlurper.parseText(file("https://zenodo.org/api/records/${zenodo_id}").text).metadata.version
    meta = [id: 'checkm2_db', version: db_version]
    """
    # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
    # So we download the database manually
    wget https://zenodo.org/records/${zenodo_id}/files/checkm2_database.tar.gz

    tar -xzf checkm2_database.tar.gz
    db_path=\$(find -name *.dmnd)
    MD5=\$(grep -o '\\.dmnd": "[^"]*"' CONTENTS.json | cut -d '"' -f 3)

    md5sum -c <<< "\$MD5  \$db_path"
    mv \$db_path checkm2_db_v${db_version}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
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
