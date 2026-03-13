process TREERECS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treerecs:1.2--h9f5acd7_3':
        'biocontainers/treerecs:1.2--h9f5acd7_3' }"

    input:
    tuple val(meta), path(species_tree), path(gene_trees)

    output:
    tuple val(meta), path("*.nwk"), emit: corrected_trees
    tuple val("${task.process}"), val("treerecs"), eval("treerecs --version 2>/dev/null | sed 's/ (.*) //g'"), topic: versions, emit: versions_treerecs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'treerecs'

    """
    mkdir -p treerecs_output

    treerecs \\
         -s ${species_tree} \\
         -g ${gene_trees} \\
         $args

    mv treerecs_output/*.nwk ${prefix}.nwk || true
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.nwk
    """
}
