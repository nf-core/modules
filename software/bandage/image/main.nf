// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BANDAGE_IMAGE {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::bandage=0.8.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bandage:0.8.1--hc9558a2_2"
    } else {
        container "quay.io/biocontainers/bandage:0.8.1--hc9558a2_2"
    }

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path('*.png'), emit: png
    tuple val(meta), path('*.svg'), emit: svg
    path  '*.version.txt'         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    Bandage image $gfa ${prefix}.png $options.args
    Bandage image $gfa ${prefix}.svg $options.args

    echo \$(Bandage --version 2>&1) | sed 's/^.*Version: //; s/ .*\$//' > ${software}.version.txt
    """
}
