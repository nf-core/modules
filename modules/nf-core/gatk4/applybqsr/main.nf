process GATK4_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1330cb37b3acde07db0a04befa76973fea72d1d10a79542132e3d77950225af/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:b945148cc5e2a616'}"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    tuple val(meta2), path(fasta), path(fai), path(dict)
    val(output_suffix)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam, optional: true
    tuple val(meta), path("${prefix}*bai"), emit: bai, optional: true
    tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // suffix can only be bam or cram, cram being the sensible default
    def suffix = output_suffix == "bam" ? "bam" : "cram"

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
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = output_suffix == "bam" ? "bam" : "cram"

    """
    touch ${prefix}.${suffix}
    if [[ ${suffix} == cram ]]; then
        touch ${prefix}.cram.bai
    else
        touch ${prefix}.bai
    fi
    """
}
