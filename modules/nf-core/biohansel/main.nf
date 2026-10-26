process BIOHANSEL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1a531965c239cfcd2a4055dbd1258e5dc24fe06da4393df52d9102de0e1100e/data'
        : 'community.wave.seqera.io/library/bio_hansel_numpy_pandas_python_setuptools:7c4db842c7382657' }"

    input:
    tuple val(meta), path(seqs)
    path scheme_metadata

    output:
    tuple val(meta), path("${prefix}-summary.txt")       , emit: summary
    tuple val(meta), path("${prefix}-kmer-results.txt")  , emit: kmer_results
    tuple val(meta), path("${prefix}-simple-summary.txt"), emit: simple_summary
    tuple val("${task.process}"), val('biohansel'), eval("hansel --version 2>&1 | sed 's/^.*hansel //'"), topic: versions, emit: versions_biohansel

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
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-summary.txt
    touch ${prefix}-kmer-results.txt
    touch ${prefix}-simple-summary.txt
    """
}
