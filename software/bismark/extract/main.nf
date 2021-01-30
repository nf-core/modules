// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BISMARK_EXTRACT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bismark==0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*.gz"), emit: consensus
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt"), emit: mbias
    path  "*.version.txt"         , emit: version

    script:
    //Assign sensible numbers for multicore and buffer_size based on bismark docs
    def ccore = task.cpus ? ((task.cpus as int) / 3) as int : 1
    def multicore = (ccore > 1) ? "--multicore ${ccore}" : ""
    //Only set buffer_size when there are more than 6.GB of memory available
    def mbuffer = task.memory ? (task.memory as nextflow.util.MemoryUnit) - 2.GB : (4.GB).toBytes()
    def buffer = (mbuffer.compareTo(4.GB) == 1) ? "--buffer_size ${mbuffer.toGiga()}G" : ""


    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def seqtype = meta.single_end ? '-s' : '-p --ignore_r2 2 --ignore_3prime_r2 2 --no_overlap'

    """
    bismark_methylation_extractor \\
        --bedGraph \\
        --counts \\
        --gzip \\
        --report \\
        $seqtype \\
        $options.args \\
        $multicore \\
        $buffer \\
        $bam

    echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//' > ${software}.version.txt
    """
}
