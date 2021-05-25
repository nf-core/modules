// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRIX {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pairix=0.3.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairix:0.3.7--py36h30a8e3e_3"
    } else {
        container "quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3"
    }

    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path(pair), path("*.px2"), emit: index
    path "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    pairix \\
        $options.args \\
        $pair

    echo \$(pairix --help 2>&1) | sed 's/^.*Version: //; s/Usage.*\$//' > ${software}.version.txt
    """
}
