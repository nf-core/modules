
process MMSEQS_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321h6a68c12_2':
        'biocontainers/mmseqs2:14.7e284--pl5321h6a68c12_2' }"

    input:
    tuple val(meta), path(db_query)
    tuple val(meta2), path(db_target)

    output:
    tuple val(meta), path("${prefix}"), emit: db_search
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if ("$db_query" == "${prefix}" || "$db_target" == "${prefix}"  ) error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p ${prefix}

    # Extract basename of DB
    DB_QUERY_PATH_NAME=\$(find -L "$db_query/" -name "*.dbtype" | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' | head -n 1 )
    DB_TARGET_PATH_NAME=\$(find -L "$db_target/" -name "*.dbtype" | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D'| head -n 1 )

    mmseqs \\
        search \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix}/${prefix} \\
        tmp1 \\
        --threads ${task.cpus} \\
        --compressed 1 \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$db_query" == "${prefix}" || "$db_target" == "${prefix}"  ) error "Input and output names of databases are the same, set prefix in module configuration to disambiguate!"
    """
    mkdir -p $prefix
    touch ${prefix}/${prefix}.{0..9}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
