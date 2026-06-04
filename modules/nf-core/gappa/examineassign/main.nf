process GAPPA_EXAMINEASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gappa:0.8.0--h9a82719_0':
        'quay.io/biocontainers/gappa:0.8.0--h9a82719_0' }"

    input:
    tuple val(meta), path(jplace), path(taxonomy)

    output:
    tuple val(meta), path("*profile.tsv")         , emit: profile
    tuple val(meta), path("*labelled_tree.newick"), emit: labelled_tree
    tuple val(meta), path("*per_query.tsv")       , emit: per_query, optional: true
    tuple val(meta), path("*krona.profile")       , emit: krona    , optional: true
    tuple val(meta), path("*sativa.tsv")          , emit: sativa   , optional: true
    tuple val("${task.process}"), val('gappa'), eval("gappa --version 2>&1 | sed 's/v//'"), emit: versions_gappa, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine assign \\
        ${args} \\
        --threads ${task.cpus} \\
        --jplace-path ${jplace} \\
        --taxon-file ${taxonomy} \\
        --file-prefix ${prefix}.
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.tsv
    touch ${prefix}.labelled_tree.newick
    """
}
