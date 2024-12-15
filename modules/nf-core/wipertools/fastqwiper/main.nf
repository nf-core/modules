process WIPERTOOLS_FASTQWIPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq_in)

    output:
    tuple val(meta), path("${prefix}.fastq.gz") , emit: fastq_out
    path("*.report")                            , emit: report
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    prefix      = prefix + "_wiped"
    """
    wipertools \\
        fastqwiper \\
        -i ${fastq_in} \\
        -o ${prefix}.fastq.gz \\
        -r ${prefix}.report \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqwiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    prefix          = prefix + "_wiped"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqwiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """
}
