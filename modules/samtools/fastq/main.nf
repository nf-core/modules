// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::samtools=1.13' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
    } else {
        container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def endedness = meta.single_end ? "-0 ${prefix}.fastq.gz" : "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz"

    """
    samtools fastq \\
        $options.args \\
        -@ $task.cpus \\
        $endedness \\
        $bam
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
