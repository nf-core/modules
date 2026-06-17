process WIPERTOOLS_FASTQGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.5--pyhdfd78af_0':
        'quay.io/biocontainers/wipertools:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.fastq.gz"), emit: gathered_fastq
    tuple val("${task.process}"), val('wipertools'), eval("wipertools fastqgather --version"), topic: versions, emit: versions_wipertools


    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_gather"
    fastq_string = fastq.collect{ file -> file.name }.sort().join(" ")

    // Check if the output file name is in the list of input files
    if (fastq.any { file -> file.name == "${prefix}.fastq.gz" }) {
        error 'Output file name "${prefix}.fastq.gz" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    wipertools \\
        fastqgather \\
        -i ${fastq_string} \\
        -o ${prefix}.fastq.gz \\
        ${args}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (fastq.any { file -> file.name == "${prefix}.fastq.gz" }) {
        error 'Output file name "${prefix}.fastq.gz" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }
    """
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
