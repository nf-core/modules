// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COOLER_ZOOMIFY {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0"
    }

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.mcool"), emit: mcool
    path "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cooler zoomify \\
        $options.args \\
        -n $task.cpus \\
        -o ${prefix}.mcool \\
        $cool

    echo \$(cooler --version 2>&1) | sed 's/cooler, version //' > ${software}.version.txt
    """
}
