process NEXTCLADE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::nextclade_js=0.14.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade_js:0.14.4--h9ee0642_0' :
        'quay.io/biocontainers/nextclade_js:0.14.4--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.csv")       , emit: csv
    tuple val(meta), path("${prefix}.json")      , emit: json
    tuple val(meta), path("${prefix}.tree.json") , emit: json_tree
    tuple val(meta), path("${prefix}.tsv")       , emit: tsv
    tuple val(meta), path("${prefix}.clades.tsv"), optional:true, emit: tsv_clades
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    nextclade \\
        $args \\
        --jobs $task.cpus \\
        --input-fasta $fasta \\
        --output-json ${prefix}.json \\
        --output-csv ${prefix}.csv \\
        --output-tsv ${prefix}.tsv \\
        --output-tsv-clades-only ${prefix}.clades.tsv \\
        --output-tree ${prefix}.tree.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1)
    END_VERSIONS
    """
}
