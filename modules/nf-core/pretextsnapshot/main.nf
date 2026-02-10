process PRETEXTSNAPSHOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextsnapshot:0.0.5--h9948957_0':
        'biocontainers/pretextsnapshot:0.0.5--h9948957_0' }"

    input:
    tuple val(meta), path(pretext_map)

    output:
    tuple val(meta), path('*.{jpeg,png,bmp}'), emit: image
    tuple val("${task.process}"), val('PretextSnapshot'), eval("PretextSnapshot --version | sed 's/^.*PretextSnapshot Version //g'"), emit: versions_pretextsnapshot, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_"
    """
    PretextSnapshot \\
        $args \\
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
