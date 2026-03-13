process SAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam")  , emit: bam ,   optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.bai")  , emit: bai ,   optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        $args \\
        $input \\
        -o ${prefix}.${output_extension}

    samtools index -@${task.cpus} ${prefix}.${output_extension}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"
    def index_extension = output_extension == "bam" ? "bai" : "crai"

    """
    touch ${prefix}.${output_extension}
    touch ${prefix}.${output_extension}.${index_extension}
    """
}
