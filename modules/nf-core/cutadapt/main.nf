process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py311haab0aaa_0' :
        'biocontainers/cutadapt:5.2--py311haab0aaa_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    tuple val("${task.process}"), val("cutadapt"), eval('cutadapt --version'), topic: versions, emit: versions_cutadapt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    """
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.trim.fastq.gz ;"
    }
    else {
        output_command  = "echo '' | gzip > ${prefix}_1.trim.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}_2.trim.fastq.gz ;"
    }
    """
    ${output_command}
    touch ${prefix}.cutadapt.log
    """
}
