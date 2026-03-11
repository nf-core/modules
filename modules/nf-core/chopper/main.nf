process CHOPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.9.0--hdcf5f25_0':
        'biocontainers/chopper:0.9.0--hdcf5f25_0' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    tuple val("${task.process}"), val('chopper'), eval("chopper --version 2>&1 | cut -d ' ' -f 2"), emit: versions_chopper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def args3  = task.ext.args3  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_filtering = fasta ? "--contam ${fasta}" : ""

    if ("$fastq" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    zcat \\
        $args \\
        $fastq | \\
    chopper \\
        --threads $task.cpus \\
        $fasta_filtering \\
        $args2 | \\
    gzip \\
        $args3 > ${prefix}.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.fastq.gz
    """
}