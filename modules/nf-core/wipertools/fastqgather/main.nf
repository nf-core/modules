process WIPERTOOLS_FASTQGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq_in)

    output:
    tuple val(meta), path("*_gather.fastq.gz") , emit: fastq_out
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    wipertools \\
        fastqgather \\
        -i $fastq_in \\
        -o ${prefix}_gather.fastq.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqgather: \$(wipertools fastqgather --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_gather.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqgather: \$(wipertools fastqgather --version)
    END_VERSIONS
    """
}
