process GATK4_SAMTOFASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b9822b92da68a3e7916072218082e3fa79bebc2f377947c363613adeecd56ec5/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:961440660027ec01'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = meta.single_end ? "--FASTQ ${prefix}.fastq.gz" : "--FASTQ ${prefix}_1.fastq.gz --SECOND_END_FASTQ ${prefix}_2.fastq.gz"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SamToFastq \\
        --INPUT ${bam} \\
        ${output} \\
        --TMP_DIR . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    echo "" | gzip > ${prefix}_1.fastq.gz
    echo "" | gzip > ${prefix}_2.fastq.gz
    """
}
