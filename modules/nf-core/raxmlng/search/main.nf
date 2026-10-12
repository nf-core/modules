process RAXMLNG_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:2.0.3--h870a6a7_0' :
        'quay.io/biocontainers/raxml-ng:2.0.3--h870a6a7_0' }"

    input:
    tuple val(meta), path(alignment), val(model)
    path(tree)
    path(tree_constraint)
    path(partitions)
    path(site_weights)

    output:
    tuple val(meta), path("*.raxml.bestTree") , emit: phylogeny
    tuple val(meta), path("*.raxml.bestModel"), emit: best_model
    tuple val(meta), path("*.raxml.mlTrees")  , emit: ml_trees  , optional: true
    tuple val(meta), path("*.raxml.startTree"), emit: start_tree, optional: true
    tuple val(meta), path("*.raxml.log")      , emit: log
    tuple val("${task.process}"), val('raxmlng'), eval("raxml-ng --version 2>&1 | sed '/RAxML-NG v/!d;s/.*v. //;s/ .*//'"), emit: versions_raxmlng, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fix random seed for reproducibility if not specified in command line
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=42"}
    def tree_arg             = tree             ? "--tree ${tree}"                       : ''
    def tree_constraint_arg  = tree_constraint  ? "--tree-constraint ${tree_constraint}"  : ''
    // a partitions/model-definition file takes the place of a literal model string;
    // --model only ever takes one or the other, never both
    def model_arg            = partitions       ? "--model ${partitions}"                : "--model ${model}"
    def site_weights_arg     = site_weights     ? "--site-weights ${site_weights}"        : ''
    """
    raxml-ng \\
        --search \\
        ${args} \\
        --msa ${alignment} \\
        ${model_arg} \\
        ${tree_arg} \\
        ${tree_constraint_arg} \\
        ${site_weights_arg} \\
        --threads ${task.cpus} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raxml.bestTree
    touch ${prefix}.raxml.bestModel
    touch ${prefix}.raxml.mlTrees
    touch ${prefix}.raxml.startTree
    touch ${prefix}.raxml.log
    """
}
