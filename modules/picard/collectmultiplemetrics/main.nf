process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf")    , emit: pdf
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectMultipleMetrics \\
        $args \\
        INPUT=$bam \\
        OUTPUT=${prefix}.CollectMultipleMetrics \\
        REFERENCE_SEQUENCE=$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectMultipleMetrics.alignment_summary_metrics
    touch ${prefix}.CollectMultipleMetrics.insert_size_metrics
    touch ${prefix}.CollectMultipleMetrics.quality_distribution.pdf
    touch ${prefix}.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    touch ${prefix}.CollectMultipleMetrics.quality_by_cycle_metrics
    touch ${prefix}.CollectMultipleMetrics.read_length_histogram.pdf
    touch ${prefix}.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    touch ${prefix}.CollectMultipleMetrics.quality_by_cycle.pdf
    touch ${prefix}.CollectMultipleMetrics.insert_size_histogram.pdf
    touch ${prefix}.CollectMultipleMetrics.quality_distribution_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
