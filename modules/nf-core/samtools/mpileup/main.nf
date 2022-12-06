process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"
    input:
    tuple val(meta), path(input), path(intervals)
    path  fasta

    output:
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = intervals ? "-l ${intervals}" : ""
    """
    samtools mpileup \\
        --fasta-ref $fasta \\
        --output ${prefix}.mpileup \\
        $args \\
        $intervals \\
        $input
    bgzip ${prefix}.mpileup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
