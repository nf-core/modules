process MMSEQS_TSV2EXPROFILEDB {
    tag "${database}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0'
        : 'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0'}"

    input:
    path database

    output:
    path (database)    , emit: db_exprofile
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    """
    DB_PATH_NAME=\$(find -L "${database}/" -name "*_seq.tsv" | sed 's/_seq\\.tsv\$//')
    touch "\${DB_PATH_NAME}_db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
