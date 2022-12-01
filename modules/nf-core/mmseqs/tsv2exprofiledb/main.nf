process MMSEQS_TSV2EXPROFILEDB {
    tag "$db"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::mmseqs2=14.7e284" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0':
        'quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0' }"

    input:
    path db
    val db_name

    output:
    path (db)          , emit: db_exprofile
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_path_name = db_name ? "${db}/${db_name}": "${db}/${db}"
    """
    mmseqs tsv2exprofiledb \\
        $db_path_name \\
        "${db_path_name}_db" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${db_path_name}_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
