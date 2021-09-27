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

    conda (params.enable_conda ? "bioconda::pbccs=6.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbccs:6.0.0--h9ee0642_2"
    } else {
        container "quay.io/biocontainers/pbccs:6.0.0--h9ee0642_2"
    }

    input:
    tuple val(meta), path(bam), path(pbi)
    val chunk_num
    val chunk_on

    output:
    tuple val(meta), path("*.ccs.bam")            , emit: bam
    tuple val(meta), path("*.ccs.bam.pbi")        , emit: pbi
    tuple val(meta), path("*.ccs_report.txt" )    , emit: ccs_report_txt
    tuple val(meta), path("*.ccs_report.json" )   , emit: ccs_report_json
    tuple val(meta), path("*.zmw_metrics.json.gz"), emit: zmw_metrics
    path "versions.yml"                           , emit: version

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ccs \\
        $bam \\
        ${prefix}.ccs.bam \\
        --report-file ${prefix}.${chunk_num}.ccs_report.txt \\
        --report-json ${prefix}.${chunk_num}.ccs_report.json \\
        --metrics-json ${prefix}.${chunk_num}.zmw_metrics.json.gz \\
        --chunk $chunk_num/$chunk_on \\
        -j $task.cpus \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ccs: \$( ccs --version|head -n 1|sed 's/ccs //'|sed 's/ (.*//' )
    END_VERSIONS
    """
}
