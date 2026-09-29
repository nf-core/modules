process GATK4_SELECTVARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1be338209b229f0bfbddf022580e1e1d9a7ebb095bc69db67c6485dc24402959/data'
        : 'community.wave.seqera.io/library/gatk4-main_gzip:2e707932c78d8e37'}"

    input:
    tuple val(meta), path(vcf), path(vcf_idx), path(intervals)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK SelectVariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --variant ${vcf} \\
        --output ${prefix}.vcf.gz \\
        ${interval} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
