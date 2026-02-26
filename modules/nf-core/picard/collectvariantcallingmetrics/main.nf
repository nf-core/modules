process PICARD_COLLECTVARIANTCALLINGMETRICS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0'
        : 'biocontainers/picard:3.4.0--hdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(index), path(intervals_file), path(fasta), path(dict), path(dbsnp), path(dbsnp_index)

    output:
    tuple val(meta), path("*.variant_calling_detail_metrics"), emit: detail_metrics
    tuple val(meta), path("*.variant_calling_summary_metrics"), emit: summary_metrics
    tuple val("${task.process}"), val('picard'), eval("picard CollectVariantCallingMetrics --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def intervals = intervals_file ? "--TARGET_INTERVALS ${intervals_file}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard CollectVariantCallingMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectVariantCallingMetrics \\
        ${args} \\
        --THREAD_COUNT ${task.cpus} \\
        --INPUT ${vcf} \\
        --OUTPUT ${prefix} \\
        --DBSNP ${dbsnp} \\
        ${reference} \\
        --TMP_DIR . \\
        ${intervals} \\
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_calling_detail_metrics
    touch ${prefix}.variant_calling_summary_metrics
    """
}
