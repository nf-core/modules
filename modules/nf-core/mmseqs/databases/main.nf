process MMSEQS_DATABASES {
    tag '$database'
    label 'process_medium'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321h6a68c12_2':
        'biocontainers/mmseqs2:14.7e284--pl5321h6a68c12_2' }"

    input:
    val database

    output:
    path "mmseqs_database/" , emit: database
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir mmseqs_database/

    mmseqs databases \\
        ${database} \\
        mmseqs_database/database \\
        tmp/ \\
        --threads ${task.cpus} \\
        --compressed 1 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir mmseqs_database/

    touch mmseqs_database/database
    touch mmseqs_database/database.dbtype
    touch mmseqs_database/database_h
    touch mmseqs_database/database_h.dbtype
    touch mmseqs_database/database_h.index
    touch mmseqs_database/database.index
    touch mmseqs_database/database.lookup
    touch mmseqs_database/database_mapping
    touch mmseqs_database/database.source
    touch mmseqs_database/database_taxonomy
    touch mmseqs_database/database.version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
