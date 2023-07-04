
process MMSEQS_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321h6a68c12_2':
        'biocontainers/mmseqs2:14.7e284--pl5321h6a68c12_2' }"

    input:
    tuple val(meta), path(query_db)
    tuple val(meta2), path(target_db)

    output:
    tuple val(meta), path("${prefix}") , emit: search_db
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$query_db" == "${prefix}" || "$target_db" == "${prefix}"  ) error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p ${prefix}

    QUERY_DB_PATH_NAME=\$(find -L "$query_db/" -name "*.dbtype" | sed 's/\\.dbtype\$//' | head -n 1 )
    TARGET_DB_PATH_NAME=\$(find -L "$target_db/" -name "*.dbtype" | sed 's/\\.dbtype\$//'| head -n 1 )

    mmseqs \\
        search \\
        \$QUERY_DB_PATH_NAME \\
        \$TARGET_DB_PATH_NAME \\
        ${prefix}/search \\
        tmp1 \\
        $args \\
        --threads ${task.cpus} \\
        --compressed 1


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$query_db" == "${prefix}" || "$target_db" == "${prefix}"  ) error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p $prefix
    touch ${prefix}/search.{0..9}
    touch ${prefix}/search.dbtype
    touch ${prefix}/search.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
