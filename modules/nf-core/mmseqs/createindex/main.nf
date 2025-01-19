process MMSEQS_CREATEINDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:16.747c6--pl5321h6a68c12_0'
        : 'biocontainers/mmseqs2:16.747c6--pl5321h6a68c12_0'}"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path(db), emit: db_indexed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        createindex \\
        \${DB_INPUT_PATH_NAME} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: "*.dbtype"
    """
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    touch "\${DB_PATH_NAME}.idx"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
