// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PBCCS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pbccs=6.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbccs:6.0.0--h9ee0642_2"
    } else {
        container "quay.io/biocontainers/pbccs:6.0.0--h9ee0642_2"
    }

    input:
    tuple val(meta), path(bam), path(pbi), val(chunk_num), val(chunk_on)

    output:
    tuple val(meta), path("*.ccs.bam"), path("*.ccs.bam.pbi"), emit: bam
    tuple path("*.ccs_report.txt"), path("*.ccs_report.json"), path ("*.zmw_metrics.json.gz"), emit: reports
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def ccs         = bam.toString().replaceAll(/bam$/, '') + chunk_num + '.ccs.bam'
    def report_txt  = bam.toString().replaceAll(/bam$/, '') + chunk_num + '.ccs_report.txt'
    def report_json = bam.toString().replaceAll(/bam$/, '') + chunk_num + '.ccs_report.json'
    def zmw_metrics = bam.toString().replaceAll(/bam$/, '') + chunk_num + '.zmw_metrics.json.gz'
    """
    ccs \\
        $bam \\
        $ccs \\
        --report-file $report_txt \\
        --report-json $report_json \\
        --metrics-json $zmw_metrics \\
        --chunk $chunk_num/$chunk_on \\
        -j $task.cpus \\
        $options.args

    echo \$(ccs --version 2>&1) | grep -e 'commit' > ${software}.version.txt
    """
}
