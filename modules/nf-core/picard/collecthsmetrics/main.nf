process PICARD_COLLECTHSMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'quay.io/biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path bait_intervals
    path target_intervals

    output:
    tuple val(meta), path("*_metrics")  , emit: metrics
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CollectHsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectHsMetrics \\
        $args \\
        $reference \\
        --BAIT_INTERVALS $bait_intervals \\
        --TARGET_INTERVALS $target_intervals \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.CollectHsMetrics.coverage_metrics


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectHsMetrics.coverage_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
