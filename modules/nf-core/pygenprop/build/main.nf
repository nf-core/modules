process PYGENPROP_BUILD {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/025a1b9aa4f0042c0af7e3a8acc54635876014f2b247ac801d91874d97440c99/data':
        'community.wave.seqera.io/library/pip_python_pygenprop:829a7d86185a9fd4' }"

    input:
    tuple val(meta), path(ips)
    path(gp_txt)

    output:
    tuple val(meta), path("*.micro"), emit: meda
    tuple val("${task.process}"), val('pygenprop'), val('1.1'), topic: versions, emit: versions_pygenprop
    tuple val("${task.process}"), val('python'), eval('python -V | sed "s/Python //g"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pygenprop \\
        build \\
        $args \\
        -d ${gp_txt} \\
        -i ${ips} \\
        -o ${prefix}.micro
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.micro
    """
}
