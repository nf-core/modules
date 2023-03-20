process VSEARCH_SORT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vsearch=2.22.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.22.1':
        'quay.io/biocontainers/vsearch:2.22.1' }"

    input:
    tuple val(meta), path(fasta)
    val sort_arg

    output:
    tuple val(meta), path("*_sorted.fasta"), emit: fasta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \
        $sort_arg $fasta \
        --output ${prefix}_sorted.fasta
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
