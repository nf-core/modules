// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BISMARK_METHYLATIONEXTRACTOR {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bismark:0.23.0--0"
    } else {
        container "quay.io/biocontainers/bismark:0.23.0--0"
    }

    input:
    tuple val(meta), path(bam)
    path index

    output:
    tuple val(meta), path("*.bedGraph.gz")         , emit: bedgraph
    tuple val(meta), path("*.txt.gz")              , emit: methylation_calls
    tuple val(meta), path("*.cov.gz")              , emit: coverage
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt")          , emit: mbias
    path "versions.yml"                            , emit: versions

    script:
    // Assign sensible numbers for multicore and buffer_size based on bismark docs
    def ccore = task.cpus ? ((task.cpus as int) / 3) as int : 1
    def multicore = (ccore > 1) ? "--multicore ${ccore}" : ""
    // Only set buffer_size when there are more than 6.GB of memory available
    def mbuffer = task.memory ? (task.memory as nextflow.util.MemoryUnit) - 2.GB : (4.GB).toBytes()
    def buffer = (mbuffer.compareTo(4.GB) == 1) ? "--buffer_size ${mbuffer.toGiga()}G" : ""

    def seqtype  = meta.single_end ? '-s' : '-p'
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

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
