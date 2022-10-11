process BIOHANSEL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bio_hansel=2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bio_hansel%3A2.6.1--py_0':
        'quay.io/biocontainers/bio_hansel:2.6.1--py_0' }"

    input:
    tuple val(meta), path(seqs)
    val scheme
    path scheme_metadata

    output:
    tuple val(meta), path("summary.txt")       , emit: summary
    tuple val(meta), path("kmer-results.txt")  , emit: kmer_results
    tuple val(meta), path("simple-summary.txt"), emit: simple_summary
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def scheme_metadata_opt = scheme_metadata ? "--scheme-metadata ${scheme_metadata[0]}" : ""
    def input_type = seqs[1] == null ? "" : "--paired-reads"
    """
    hansel \\
        $args \\
        $scheme_metadata_opt \\
        --threads $task.cpus \\
        --scheme $scheme \\
        --output-summary summary.txt \\
        --output-kmer-results kmer-results.txt \\
        --output-simple-summary simple-summary.txt \\
        $input_type $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biohansel: \$(echo \$(hansel --version 2>&1) | sed 's/^.*hansel //' )
    END_VERSIONS
    """
}
