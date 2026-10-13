process GATK4_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b33544fcc9268dad01801b928e566f1beaacd078b1d58d7904d7c6aa68eba9a/data'
        : 'community.wave.seqera.io/library/gatk4-lite_samtools:9b6806e26dd32afe'}"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    tuple val(meta2), path(fasta), path(fai), path(dict)
    val output_suffix

    output:
    tuple val(meta), path("*cram"), emit: cram, optional: true
    tuple val(meta), path("*bam"), emit: bam, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // suffix can only be bam or cram, cram being the sensible default
    def suffix = output_suffix == "bam" ? "bam" : "cram"
    def gatk_index = suffix == "cram" ? "-OBI false" : ""

    def interval_command = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input ${input} \\
        --output ${prefix}.${suffix} \\
        --reference ${fasta} \\
        --bqsr-recal-file ${bqsr_table} \\
        ${interval_command} \\
        ${gatk_index} \\
        --tmp-dir . \\
        ${args}

    if [[ ${suffix} == cram ]]; then
        samtools index ${prefix}.cram
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = output_suffix == "bam" ? "bam" : "cram"

    """
    touch ${prefix}.${suffix}
    if [[ ${suffix} == cram ]]; then
        touch ${prefix}.cram.crai
    else
        touch ${prefix}.bai
    fi
    """
}
