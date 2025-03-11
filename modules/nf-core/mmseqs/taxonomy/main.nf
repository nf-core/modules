process MMSEQS_TAXONOMY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    tuple val(meta), path(db_query)
    path(db_target)

    output:
    tuple val(meta), path("${prefix}_taxonomy"), emit: db_taxonomy
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype" //represents the db_query
    def args3 = task.ext.args3 ?: "*.dbtype" //represents the db_target
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_taxonomy

    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_QUERY_PATH_NAME=\$(find -L "${db_query}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_TARGET_PATH_NAME=\$(find -L "${db_target}/" -maxdepth 1 -name "${args3}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        taxonomy \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix}_taxonomy/${prefix} \\
        tmp1 \\
        $args \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_taxonomy
    touch ${prefix}_taxonomy/${prefix}.{0..25}
    touch ${prefix}_taxonomy/${prefix}.dbtype
    touch ${prefix}_taxonomy/${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
