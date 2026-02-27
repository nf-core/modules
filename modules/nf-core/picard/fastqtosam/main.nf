process PICARD_FASTQTOSAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('picard'), eval("picard FastqToSam --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        log.warn('[Picard FastqToSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    def input = meta.single_end ? "--FASTQ ${reads}" : "--FASTQ ${reads[0]} --FASTQ2 ${reads[1]}"
    def sample_name = args.contains("--SAMPLE_NAME") || args.contains("-SM") ? "" : "--SAMPLE_NAME ${prefix}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        FastqToSam \\
        ${args} \\
        ${input} \\
        ${sample_name} \\
        --OUTPUT ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
