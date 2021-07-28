// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PICARD_MERGESAMFILES {
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
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bam_files = bams.sort()
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    if (bam_files.size() > 1) {
        """
        picard \\
            -Xmx${avail_mem}g \\
            MergeSamFiles \\
            $options.args \\
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${prefix}.bam
        echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
        """
    } else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        echo \$(picard MergeSamFiles --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d: > ${software}.version.txt
        """
    }
}
