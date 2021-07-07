// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BISMARK_DEDUPLICATE {
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

    output:
    tuple val(meta), path("*.deduplicated.bam")        , emit: bam
    tuple val(meta), path("*.deduplication_report.txt"), emit: report
    path  "*.version.txt"                              , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def seqtype    = meta.single_end ? '-s' : '-p'
    """
    deduplicate_bismark \\
        $options.args \\
        $seqtype \\
        --bam $bam

    echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//' > ${software}.version.txt
    """
}
