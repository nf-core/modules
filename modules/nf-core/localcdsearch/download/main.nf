process LOCALCDSEARCH_DOWNLOAD {
    tag "${databases.join(', ')}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/local-cd-search:0.3.0--pyhdfd78af_0' :
        'biocontainers/local-cd-search:0.3.0--pyhdfd78af_0' }"

    input:
    path download_path
    val databases

    output:
    path download_path, emit: db
    tuple val("${task.process}"), val('local-cd-search'), eval("echo '0.3.0'"), topic: versions, emit: versions_localcdsearch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    local-cd-search \\
        download \\
        ${args} \\
        ${download_path} \\
        ${databases.join(' ')}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    """
}
