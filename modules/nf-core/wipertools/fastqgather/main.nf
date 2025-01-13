process WIPERTOOLS_FASTQGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)    // channel: [ val(meta), [ .fastq|.fastq.gz ] ]

    output:
    tuple val(meta), path("${prefix}.fastq.gz"), emit: gathered_fastq   // channel: [ val(meta), *_gather.fastq | *_gather.fastq.gz ]
    path "versions.yml"                        , emit: versions         // channel: [ versions.yml ]

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (fastq.any { it.name == "${prefix}.fastq.gz" }) {
        error 'Output file name "${prefix}.fastq.gz}" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    wipertools \\
        fastqgather \\
        -i $fastq \\
        -o ${prefix}.fastq.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqgather: \$(wipertools fastqgather --version)
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (fastq.any { it.name == "${prefix}.fastq.gz" }) {
        error 'Output file name "${prefix}.fastq.gz}" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }
    """
    echo "" | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqgather: \$(wipertools fastqgather --version)
    END_VERSIONS
    """
}
