process RESEQ_ILLUMINAPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/reseq:1.1--py310hfb68e69_5' :
        'biocontainers/reseq:1.1--py310hfb68e69_5' }"

    input:
    tuple val(meta), path(bam), path(fasta), path(adapter_fasta), path(adapter_matrix)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fastq
    tuple val("${task.process}"), val('reseq'), eval('reseq --version 2>&1 | sed \'s/^.*ReSeq version //\''), emit: versions_reseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_arg = adapter_fasta ? "--adapterFile ${adapter_fasta}" : ''
    def matrix_arg = adapter_matrix ? "--adapterMatrix ${adapter_matrix}" : ''
    """
    reseq \\
        illuminaPE \\
        -j $task.cpus \\
        -r $fasta \\
        -b $bam \\
        $adapter_arg \\
        $matrix_arg \\
        $args \\
        -1 ${prefix}_R1.fq \\
        -2 ${prefix}_R2.fq

    gzip ${prefix}_R1.fq
    gzip ${prefix}_R2.fq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}_R1.fq.gz
    echo '' | gzip > ${prefix}_R2.fq.gz
    """
}
