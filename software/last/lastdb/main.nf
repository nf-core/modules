// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_LASTDB {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::last=1219" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1219--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1219--h2e03b76_0"
    }

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("lastdb"), emit: index
    path "*.version.txt"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    
    """
    mkdir lastdb
    lastdb \\
        ${options.args} \\
        -P ${task.cpus} \\
        lastdb/${meta.id} \\
        ${fastx}

    echo \$(lastdb --version 2>&1) | sed 's/lastdb //' > ${software}.version.txt
    """
}
