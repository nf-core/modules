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
    tuple val(meta), path("*_wiped.fastq.gz") , emit: fastq_out
    path("*.report")                          , emit: report
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def fastq_out   = prefix.endsWith('.fastq.gz') ? prefix.replaceAll(/\.fastq.gz$/, '_wiped.fastq.gz') : (prefix.endsWith('.fastq') ? prefix.replaceAll(/\.fastq$/, '_wiped.fastq.gz') : prefix + "_wiped.fastq.gz")
    def report_file = prefix + ".report"
    """
    wipertools \\
        fastqwiper \\
        -i ${fastq_in} \\
        -o ${fastq_out} \\
        -r ${report_file} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools wiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_wiped.fastq.gz
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools wiper: \$(wipertools fastqwiper --version)
    END_VERSIONS
    """
}
