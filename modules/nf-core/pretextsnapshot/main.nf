process PRETEXTSNAPSHOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fde9c7892c9878e4aa8cbcb744039045a182efabac3e311a7fdb1e43b1fbb4d6/data':
        'community.wave.seqera.io/library/pretextsnapshot:0.0.7--9470b2ea6b8991c8' }"

    input:
    tuple val(meta), path(pretext_map), path(order_file)

    output:
    tuple val(meta), path('*.{jpeg,png,bmp}'), emit: image
    tuple val("${task.process}"), val('PretextSnapshot'), eval("PretextSnapshot --version | sed 's/^.*PretextSnapshot Version //g'"), emit: versions_pretextsnapshot, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_"
    def order_arg = order_file ? "--order ${order_file}" : ""
    """
    PretextSnapshot \\
        ${args} \\
        ${order_arg} \\
        --map $pretext_map \\
        --prefix $prefix \\
        --folder .
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_"
    """
    touch ${prefix}scaffold_{1,2,3,4}.png
    """
}
