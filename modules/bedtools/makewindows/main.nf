// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_MAKEWINDOWS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1"
    }

    input:
    tuple val(meta), path(regions)
    val(use_bed)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def arg_input = use_bed ? "-b $regions" : "-g $regions"
    """
    bedtools \\
        makewindows \\
        ${arg_input} \\
        $options.args \\
        > ${prefix}.tab

    echo \$(bedtools --version) | sed -e "s/bedtools v//g" > ${software}.version.txt
    """
}
