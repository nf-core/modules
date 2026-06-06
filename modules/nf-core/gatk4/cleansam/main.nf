process GATK4_CLEANSAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_index)
    // input file must be sorted for index to be created
    val create_index
    val create_md5

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.md5"), emit: md5, optional: true
    tuple val("${task.process}"), val('gatk'), eval("gatk CleanSam --version | grep GATK | sed 's/.*(GATK) v//'"), topic: versions, emit: versions_gatk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index = create_index ? "--CREATE_INDEX true" : ""
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def md5 = create_md5 ? "--CREATE_MD5_FILE true" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CleanSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData"\\
        CleanSam \\
        ${args} \\
        ${reference} \\
        ${index} \\
        ${md5} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index = create_index ? "touch ${prefix}.bam.bai" : ""
    def md5 = create_md5 ? "touch ${prefix}.md5" : ""
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam \\
    ${index} \\
    ${md5}
    """
}
