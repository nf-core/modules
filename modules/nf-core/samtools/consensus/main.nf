process SAMTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.fasta") , emit: fasta , optional: true
    tuple val(meta), path("*.fastq") , emit: fastq , optional: true
    tuple val(meta), path("*.pileup"), emit: pileup, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-f fastq") ? "fastq" :
                    args.contains("-f pileup") ? "pileup" :
                    args.contains("-f fasta") ? "fasta" :
                    "fasta"

    """
    samtools \\
        consensus \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.${extension} \\
        $input
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-f fastq") ? "fastq" :
                    args.contains("-f pileup") ? "pileup" :
                    args.contains("-f fasta") ? "fasta" :
                    "fasta"

    """
    touch ${prefix}.${extension}
    """
}
