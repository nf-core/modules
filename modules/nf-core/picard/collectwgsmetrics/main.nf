process PICARD_COLLECTWGSMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::picard=2.27.4 r::r-base"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    path  intervallist

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    def interval  = intervallist ? "--INTERVALS ${intervallist}" : ''
    if (!task.memory) {
        log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectWgsMetrics \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.CollectWgsMetrics.coverage_metrics \\
        --REFERENCE_SEQUENCE ${fasta} \\
        $interval


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectWgsMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CollectWgsMetrics.coverage_metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectWgsMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
