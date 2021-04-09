// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALEVINQC {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bioconductor-alevinqc=1.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-alevinqc:1.6.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-alevinqc:1.6.1--r40hdfd78af_0"
    }

    input:
    tuple val(meta), path(alevin_results)

    output:
    tuple val(meta), path("alevinReport.html"), emit: report
    path "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    require(alevinQC)

    alevinQCReport(baseDir = "${alevin_results}", sampleId = "${prefix}",
                outputFile = "alevinReport.html",
                outputFormat = "html_document",
                outputDir = "./", forceOverwrite = TRUE)

    write(as.character(packageVersion("alevinQC")), paste0("${software}", ".version.txt"))
    """
}
