process WIPERTOOLS_FASTQWIPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.4--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.fastq.gz"), emit: wiped_fastq
    tuple val(meta), path("*.report")          , emit: report
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_wiped"
    if ("${fastq}" == "${prefix}.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!."
    """
    wipertools \\
        fastqwiper \\
        -i ${fastq} \\
        -o ${prefix}.fastq.gz \\
        -r ${prefix}.report \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqwiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_wiped"
    if ("${fastq}" == "${prefix}.fastq.gz") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!."
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools fastqwiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """
}
