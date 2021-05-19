// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LAST_TRAIN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::last=1219" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/last:1219--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/last:1219--h2e03b76_0"
    }

    input:
    tuple val(meta), path(fastx)
    path  index

    output:
    tuple val(meta), path("*.par"), emit: param_file
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    INDEX_NAME=`find -L lastdb/ -name "*.bck" | sed 's/.bck//' | sed 's,lastdb/,,'`

    last-train \\
        $options.args \\
        -P $task.cpus \\
        ${index}/\$INDEX_NAME \\
        $fastx \\
        > ${prefix}.\$INDEX_NAME.par

    lastdb --version | sed 's/lastdb //' > ${software}.version.txt
    """
}
