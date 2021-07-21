// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PICARD_FILTERSAMREADS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::picard=2.25.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picard:2.25.6--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/picard:2.25.6--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)
    val filter
    path readlist

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard FilterSamReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    if ( filter == 'includeAligned' || filter == 'excludeAligned' ) {
        """
        picard \\
            FilterSamReads \\
            -Xmx${avail_mem}g \\
            --INPUT $bam \\
            --OUTPUT ${prefix}.bam \\
            --FILTER $filter \\
            $options.args

        echo \$(picard FilterSamReads --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
        """
    } else if ( filter == 'includeReadList' || filter == 'excludeReadList' ) {
        """
        picard \\
            FilterSamReads \\
            -Xmx${avail_mem}g \\
            --INPUT $bam \\
            --OUTPUT ${prefix}.bam \\
            --FILTER $filter \\
            --READ_LIST_FILE $readlist \\
            $options.args

        echo \$(picard FilterSamReads --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
        """
    }
}
