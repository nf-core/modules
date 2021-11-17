// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PBCCS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pbccs=6.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbccs:6.2.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/pbccs:6.2.0--h9ee0642_0"
    }

    input:
    tuple val(meta), path(bam), path(pbi)
    val chunk_num
    val chunk_on

    output:
    tuple val(meta), path("*.chunk*.bam")     , emit: bam
    tuple val(meta), path("*.chunk*.bam.pbi") , emit: pbi
    tuple val(meta), path("*.report.txt" )    , emit: report_txt
    tuple val(meta), path("*.report.json" )   , emit: report_json
    tuple val(meta), path("*.metrics.json.gz"), emit: metrics
    path  "versions.yml"                      , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ccs \\
        $bam \\
        ${prefix}.chunk${chunk_num}.bam \\
        --report-file ${prefix}.chunk${chunk_num}.report.txt \\
        --report-json ${prefix}.chunk${chunk_num}.report.json \\
        --metrics-json ${prefix}.chunk${chunk_num}.metrics.json.gz \\
        --chunk $chunk_num/$chunk_on \\
        -j $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(ccs --version 2>&1) | grep 'ccs' | sed 's/^.*ccs //; s/ .*\$//')
    END_VERSIONS
    """
}
