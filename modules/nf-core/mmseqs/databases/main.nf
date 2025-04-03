process MMSEQS_DATABASES {
    tag "${database}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    val database

    output:
    path "${prefix}/"   , emit: database
    path "versions.yml" , emit: versions

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
