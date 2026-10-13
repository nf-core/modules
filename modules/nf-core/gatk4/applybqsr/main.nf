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
    val output_format

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: reads
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = output_format[0]
    def index_fmt = output_format[1]

    def valid_formats = [
        'bam'  : ['bai', 'csi'],
        'cram' : ['crai']
    ]
    if (!valid_formats.containsKey(file_type)) {
        error("GATK4_APPLYBQSR: invalid file type '${file_type}'. Must be 'bam' or 'cram'.")
    }
    if (!valid_formats[file_type].contains(index_fmt)) {
        error("GATK4_APPLYBQSR: invalid index format '${index_fmt}' for '${file_type}' output. Must be one of: ${valid_formats[file_type]}.")
    }

    def samtools_index_args = index_fmt == 'csi' ? '-c' : ''

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
        --output ${prefix}.${file_type} \\
        --reference ${fasta} \\
        --bqsr-recal-file ${bqsr_table} \\
        ${interval_command} \\
        -OBI false \\
        --tmp-dir . \\
        ${args}

    samtools index ${samtools_index_args} ${prefix}.${file_type}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = output_format[0]
    def index_fmt = output_format[1]

    def valid_formats = [
        'bam'  : ['bai', 'csi'],
        'cram' : ['crai']
    ]
    if (!valid_formats.containsKey(file_type)) {
        error("GATK4_APPLYBQSR: invalid file type '${file_type}'. Must be 'bam' or 'cram'.")
    }
    if (!valid_formats[file_type].contains(index_fmt)) {
        error("GATK4_APPLYBQSR: invalid index format '${index_fmt}' for '${file_type}' output. Must be one of: ${valid_formats[file_type]}.")
    }

    """
    touch ${prefix}.${file_type}
    touch ${prefix}.${file_type}.${index_fmt}
    """
}
