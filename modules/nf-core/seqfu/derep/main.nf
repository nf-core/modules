process SEQFU_DEREP {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.20.3--h1eb128b_2':
        'biocontainers/seqfu:1.20.3--h1eb128b_2' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*_derep.fasta.gz"), emit: fasta
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_derep"
    def fasta_files = fastas.collect { it.getName() }
    if (fasta_files.any { it == "${prefix}.fasta.gz" }) {
        error "Input file name coincides with the output file name: ${prefix}.fasta.gz. Please set a unique prefix."
    }

    """
    seqfu \\
        derep \\
        $args \\
        $fastas | gzip -c > "${prefix}.fasta.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_derep"
    """
    echo "" | gzip -c > "${prefix}.fasta.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """
}
