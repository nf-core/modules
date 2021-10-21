// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENRICH {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::genrich=0.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/genrich:0.6.1--h5bf99c6_1"
    } else {
        container "quay.io/biocontainers/genrich:0.6.1--h5bf99c6_1"
    }

    input:
    tuple val(meta), path(treatment_bam)
    path  control_bam
    path  blacklist_bed

    output:
    tuple val(meta), path("*narrowPeak")                     , emit: peaks
    tuple val(meta), path("*pvalues.bedGraph"), optional:true, emit: bedgraph_pvalues
    tuple val(meta), path("*pileup.bedGraph") , optional:true, emit: bedgraph_pileup
    tuple val(meta), path("*intervals.bed")   , optional:true, emit: bed_intervals
    tuple val(meta), path("*duplicates.txt")  , optional:true, emit: duplicates
    path "versions.yml"                                      , emit: versions

    script:
    def prefix     = options.suffix              ? "${meta.id}${options.suffix}"   : "${meta.id}"
    def control    = params.control_bam          ? "-c $control_bam"               : ''
    def pvalues    = params.pvalues              ? "-f ${prefix}.pvalues.bedGraph" : ""
    def pileup     = params.pileup               ? "-k ${prefix}.pileup.bedGraph"  : ""
    def bed        = params.bed                  ? "-b ${prefix}.intervals.bed"    : ""
    def duplicates = options.args.contains('-r') ? "-R ${prefix}.duplicates.txt"   : ''
    def blacklist  = params.blacklist_bed        ? "-E $blacklist_bed"             : ""
    """
    Genrich \\
        -t $treatment_bam \\
        $options.args \\
        $control \\
        $blacklist \\
        -o ${prefix}.narrowPeak \\
        $pvalues \\
        $pileup \\
        $bed \\
        $duplicates \\
        $blacklist \\
        $control

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(Genrich --version 2>&1) | sed 's/^Genrich, version //; s/ .*\$//')
    END_VERSIONS
    """
}
