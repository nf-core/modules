// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS_AMPLICONCLIP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
    } else {
        container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    }

    input:
    tuple val(meta), path(bam)
    path bed
    val save_cliprejects
    val save_clipstats

    output:
    tuple val(meta), path("*.bam")            , emit: bam
    tuple val(meta), path("*.clipstats.txt")  , optional:true, emit: stats
    tuple val(meta), path("*.cliprejects.bam"), optional:true, emit: rejects_bam
    path "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def rejects  = save_cliprejects ? "--rejects-file ${prefix}.cliprejects.bam" : ""
    def stats    = save_clipstats   ? "-f ${prefix}.clipstats.txt"               : ""
    """
    samtools \\
        ampliconclip \\
        $options.args \\
        -@ $task.cpus \\
        $rejects \\
        $stats \\
        -b $bed \\
        -o ${prefix}.bam \\
        $bam

    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
