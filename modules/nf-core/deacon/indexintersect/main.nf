process DEACON_INDEXINTERSECT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.15.0--hdd79491_0':
        'quay.io/biocontainers/deacon:0.15.0--hdd79491_0' }"

    input:
    tuple val(meta), path(indices)  // two or more deacon .idx index files to intersect

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.intersect"

    """
    deacon \\
        index \\
        intersect \\
        ${args} \\
        ${indices.join(' ')} \\
        > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    """
}
