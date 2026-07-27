process FGUMI_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a62b457c53300603da026225f95b4db04d1c9f8ba7f734787818fc105d51323/data':
        'community.wave.seqera.io/library/fgumi:0.4.0--1fb5dc6de05ce63b' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi \\
        fastq \\
        --input ${bam} \\
        --threads ${task.cpus} \\
        ${args} \\
        | gzip -c ${args2} > ${prefix}.fastq.gz
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
