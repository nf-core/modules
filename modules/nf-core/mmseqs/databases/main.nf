process MMSEQS_DATABASES {
    tag "${database}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    val database

    output:
    path "${prefix}/", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir ${prefix}/

    mmseqs databases \\
        ${database} \\
        ${prefix}/database \\
        tmp/ \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir ${prefix}/

    touch ${prefix}/database
    touch ${prefix}/database.dbtype
    touch ${prefix}/database_h
    touch ${prefix}/database_h.dbtype
    touch ${prefix}/database_h.index
    touch ${prefix}/database.index
    touch ${prefix}/database.lookup
    touch ${prefix}/database_mapping
    touch ${prefix}/database.source
    touch ${prefix}/database_taxonomy
    touch ${prefix}/database.version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
