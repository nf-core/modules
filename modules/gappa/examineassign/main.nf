process GAPPA_EXAMINEASSIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gappa=0.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gappa:0.8.0--h9a82719_0':
        'quay.io/biocontainers/gappa:0.8.0--h9a82719_0' }"

    input:
    tuple val(meta), path(jplace)
    path  taxonomy

    output:
    tuple val(meta), path("./.")                  , emit: examineassign
    tuple val(meta), path("*profile.tsv")         , emit: profile
    tuple val(meta), path("*labelled_tree.newick"), emit: labelled_tree
    tuple val(meta), path("*per_query.tsv")       , emit: per_query, optional: true
    path "versions.yml"                           , emit: versions

    when:
    taxonomy && ( task.ext.when == null || task.ext.when )

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine assign \\
        $args \\
        --threads $task.cpus \\
        --jplace-path $jplace \\
        --taxon-file $taxonomy \\
        --file-prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gappa: \$(echo \$(gappa --version 2>&1 | sed 's/v//' ))
    END_VERSIONS
    """
}
