// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::picard=2.23.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picard:2.23.9--0"
    } else {
        container "quay.io/biocontainers/picard:2.23.9--0"
    }

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf")    , emit: pdf
    path  "*.version.txt"             , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
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
        $options.args \\
        INPUT=$bam \\
        OUTPUT=${prefix}.CollectMultipleMetrics \\
        REFERENCE_SEQUENCE=$fasta

    echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version.*' | cut -f2- -d: > ${software}.version.txt
    """
}
