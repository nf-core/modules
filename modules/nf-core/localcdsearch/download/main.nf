process LOCALCDSEARCH_DOWNLOAD {
    tag "${databases.join(', ')}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/local-cd-search:0.3.0--pyhdfd78af_0' :
        'biocontainers/local-cd-search:0.3.0--pyhdfd78af_0' }"

    input:
    val databases

    output:
    path('database/'), emit: db
    tuple val("${task.process}"), val('local-cd-search'), eval("echo ${VERSION}"), topic: versions, emit: versions_localcdsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir database/
    local-cd-search \\
        download \\
        ${args} \\
        database/ \\
        ${databases.join(' ')}
    """

    stub:
    def args = task.ext.args ?: ''
    VERSION = '0.3.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo $args
    mkdir database/
    """
}
