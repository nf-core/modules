process BACPHLIP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bacphlip=0.9.6 conda-forge::numpy=1.23.5 bioconda::hmmer=3.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0':
        'biocontainers/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bacphlip")         , emit: bacphlip_results
    tuple val(meta), path("*.hmmsearch.tsv")    , emit: hmmsearch_results
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.9.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bacphlip \\
        -i $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.9.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${fasta}.bacphlip
    touch ${fasta}.hmmsearch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: $VERSION
    END_VERSIONS
    """
}
