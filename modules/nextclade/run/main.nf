process NEXTCLADE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::nextclade=1.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:1.9.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:1.9.0--h9ee0642_0' }"

    input:
    path(dataset)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.csv")       , emit: csv
    tuple val(meta), path("${prefix}.json")      , emit: json
    tuple val(meta), path("${prefix}.tree.json") , emit: json_tree
    tuple val(meta), path("${prefix}.tsv")       , emit: tsv
    tuple val(meta), path("${prefix}.clades.tsv"), optional:true, emit: tsv_clades
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    nextclade \\
        $args \\
        --jobs $task.cpus \\
        --input-fasta $fasta \\
        --input-dataset $dataset \\
        --output-json ${prefix}.json \\
        --output-csv ${prefix}.csv \\
        --output-tsv ${prefix}.tsv \\
        --output-tree ${prefix}.tree.json \\
        --output-basename ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1)
    END_VERSIONS
    """
}
