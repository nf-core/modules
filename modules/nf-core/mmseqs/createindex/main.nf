process MMSEQS_CREATEINDEX {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

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
    """
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        createindex \\
        \${DB_INPUT_PATH_NAME} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus} \\
        --split-memory-limit ${(task.memory.toGiga() * 0.8) as int}G

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: "*.dbtype"
    """
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    touch "\${DB_INPUT_PATH_NAME}.idx"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
