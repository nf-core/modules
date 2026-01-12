process LOCALCDSEARCH_DOWNLOAD {
    tag "${databases.join(', ')}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1ed921c933d8eeeb0db6d72ece09ec25edab9ad441c84b070acff1592af2d54/data' :
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
    """
    echo $args
    mkdir database/
    """
}
