// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADAM_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::adam=0.35.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/adam:0.35.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/adam:0.35.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    adam-submit \\
        --master local[${task.cpus}] \\
        --driver-memory ${task.memory.toGiga()}g \\
        -- \\
        transformAlignments \\
        -mark_duplicate_reads \\
        -single \\
        -stringency LENIENT \\
        $bam \\
        ${prefix}.md.bam

    echo \$(adam-submit --version 2>&1) | grep -o 'ADAM version: .*' | cut -f2 -d ' ' > ${software}.version.txt
    """
}
