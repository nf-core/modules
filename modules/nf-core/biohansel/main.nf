process BIOHANSEL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bio_hansel:2.6.1--py_0':
        'biocontainers/bio_hansel:2.6.1--py_0' }"

    input:
    tuple val(meta), path(seqs)
    path scheme_metadata

    output:
    tuple val(meta), path("${prefix}-summary.txt")       , emit: summary
    tuple val(meta), path("${prefix}-kmer-results.txt")  , emit: kmer_results
    tuple val(meta), path("${prefix}-simple-summary.txt"), emit: simple_summary
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def scheme_metadata_opt = scheme_metadata ? "--scheme-metadata ${scheme_metadata}" : ""
    def input_type = seqs[1] == null ? "" : "--paired-reads"
    """
    hansel \\
        $args \\
        $scheme_metadata_opt \\
        --threads $task.cpus \\
        --output-summary ${prefix}-summary.txt \\
        --output-kmer-results ${prefix}-kmer-results.txt \\
        --output-simple-summary ${prefix}-simple-summary.txt \\
        $input_type \\
        $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biohansel: \$(echo \$(hansel --version 2>&1) | sed 's/^.*hansel //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-summary.txt
    touch ${prefix}-kmer-results.txt
    touch ${prefix}-simple-summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biohansel: \$(echo \$(hansel --version 2>&1) | sed 's/^.*hansel //' )
    END_VERSIONS
    """
}
