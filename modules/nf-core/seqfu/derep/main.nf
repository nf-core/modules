process SEQFU_DEREP {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqfu:1.20.3--h1eb128b_2'
        : 'biocontainers/seqfu:1.20.3--h1eb128b_2'}"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*_derep.fasta.gz"), emit: fasta
    tuple val("${task.process}"), val('seqfu'), eval('seqfu version'), emit: versions_seqfu, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_derep"
    def fasta_files = fastas.collect { fasta_file -> fasta_file.getName() }
    if (fasta_files.any { fasta_file -> fasta_file == "${prefix}.fasta.gz" }) {
        error("Input file name coincides with the output file name: ${prefix}.fasta.gz. Please set a unique prefix.")
    }

    """
    seqfu \\
        derep \\
        ${args} \\
        ${fastas} | gzip -c > "${prefix}.fasta.gz"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_derep"
    """
    echo "" | gzip -c > "${prefix}.fasta.gz"
    """
}
