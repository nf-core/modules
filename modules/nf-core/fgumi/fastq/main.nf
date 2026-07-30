process FGUMI_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/9293544e5bd961eb471b82d28a66cdfca73fe7e70eaa504a4ce59a1b1c6e912d/data':
        'community.wave.seqera.io/library/fgumi_gzip:b3703fa2b3e8e632' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | sed "s/^fgumi //"'), topic: versions, emit: versions_fgumi
    tuple val("${task.process}"), val('gzip'), eval("gzip --version 2>&1 | sed '1!d;s/gzip //'"), emit: versions_gzip, topic: versions

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
    echo ${args}
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
