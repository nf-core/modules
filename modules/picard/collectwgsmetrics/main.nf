process PICARD_COLLECTWGSMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
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
        INPUT=$bam \\
        OUTPUT=${prefix}.CollectWgsMetrics.coverage_metrics \\
        REFERENCE_SEQUENCE=$fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectWgsMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
