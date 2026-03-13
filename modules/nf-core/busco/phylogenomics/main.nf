process BUSCO_PHYLOGENOMICS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco_phylogenomics:20240919--pyhdfd78af_0':
        'biocontainers/busco_phylogenomics:20240919--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(busco)

    output:
    tuple val(meta), path("${prefix}/gene_trees_single_copy/"), emit: gene_trees
    tuple val(meta), path("${prefix}/supermatrix/")           , emit: supermatrix
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('busco_phylogenomics'), eval('echo 20240919'), emit: versions_buscophylogenomics, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    BUSCO_phylogenomics.py \\
        -i ${busco} \\
        -o ${prefix} \\
        -t $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    mkdir ${prefix}/gene_trees_single_copy
    mkdir ${prefix}/supermatrix

    touch ${prefix}/gene_trees_single_copy/ALL.tree
    touch ${prefix}/supermatrix/SUPERMATRIX.fasta
    touch ${prefix}/supermatrix/SUPERMATRIX.partitions.nex
    touch ${prefix}/supermatrix/SUPERMATRIX.phylip
    """
}
