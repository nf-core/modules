process BACPHLIP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bacphlip=0.9.6"
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
    """
    bacphlip \\
        -i $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: 0.9.6
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${fasta}.bacphlip
    touch ${fasta}.hmmsearch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: 0.9.6
    END_VERSIONS
    """
}
