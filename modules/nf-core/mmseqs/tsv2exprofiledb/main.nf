process MMSEQS_TSV2EXPROFILEDB {
    tag "$db"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mmseqs2=14.7e284" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    path db

    output:
    path (db)          , emit: db_exprofile
    path "versions.yml", emit: versions
    // INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    DB_PATH_NAME=\$(find -L "$db/" -name "*_seq.tsv" | sed 's/_seq\\.tsv\$//')

    mmseqs tsv2exprofiledb \\
        \${DB_PATH_NAME} \\
        "\${DB_PATH_NAME}_db" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    """
    DB_PATH_NAME=\$(find -L "$db/" -name "*_seq.tsv" | sed 's/_seq\\.tsv\$//')

    touch "\${DB_PATH_NAME}_db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
