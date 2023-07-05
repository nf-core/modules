
process MMSEQS_CREATETSV {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mmseqs2=14.7e284"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(db_result)
    tuple val(meta2), path(db_query)
    tuple val(meta3), path(db_target)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    DB_RESULT_PATH_NAME=\$(find -L "$db_result/" -maxdepth 1 -name "*.dbtype" | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_QUERY_PATH_NAME=\$(find -L "$db_query/" -maxdepth 1 -name "*.dbtype" | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_TARGET_PATH_NAME=\$(find -L "$db_target/" -maxdepth 1 -name "*.dbtype" | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        createtsv \\
        ${DB_QUERY_PATH_NAME} \\
        ${DB_TARGET_PATH_NAME} \\
        ${DB_RESULT_PATH_NAME} \\
        ${prefix}.tsv \\
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
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
