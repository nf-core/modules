process TREERECS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treerecs:1.2--h9f5acd7_3' :
        'biocontainers/treerecs:1.2--h9f5acd7_3' }"

    input:
    tuple val(meta), path(gene_trees)
    path(species_tree)
    path(smap)

    output:
    tuple val(meta), path("*.nwk"), emit: corrected_trees_newick
    tuple val(meta), path("*.nhx"), emit: corrected_trees_nhx, optional: true
    tuple val(meta), path("*.phylo.xml"), emit: corrected_trees_phyloxml, optional: true
    tuple val(meta), path("*.recphylo.xml"), emit: corrected_trees_recphyloxml, optional: true
    tuple val(meta), path("*.svg"), emit: corrected_trees_svg, optional: true
    tuple val(meta), path("*.relationships_summary.txt"), emit: relationships_summary, optional: true
    tuple val("${task.process}"), val("treerecs"), eval("treerecs --version | sed -n 's/.*(\\([0-9.]*\\)).*/\\1/p'"), emit: versions_treerecs, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta?.id ?: "output"
    def smap_arg = smap ? "-S ${smap}" : ''

    """
    treerecs \\
        -s ${species_tree} \\
        -g ${gene_trees} \\
        -o treerecs_output/ \\
        ${smap_arg} \\
        ${args}

    mv treerecs_output/* . || true
    """

    stub:
    def prefix = task.ext.prefix ?: meta?.id ?: "stub"
    """
    touch ${prefix}.nwk
    touch ${prefix}.nhx
    touch ${prefix}.phylo.xml
    touch ${prefix}.recphylo.xml
    touch ${prefix}.svg
    touch ${prefix}.relationships_summary.txt

    mv treerecs_output/* . || true
    """
}