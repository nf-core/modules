// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRTOOLS_DEDUP {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pairtools=0.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5"
    } else {
        container "quay.io/biocontainers/pairtools:0.3.0--py37hb9c2fc3_5"
    }

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pairs.gz")  , emit: pairs
    tuple val(meta), path("*.pairs.stat"), emit: stat
    path "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pairtools dedup \\
        $options.args \\
        -o ${prefix}.pairs.gz \\
        --output-stats ${prefix}.pairs.stat \\
        $input

    echo \$(pairtools --version 2>&1) | sed 's/pairtools.*version //' > ${software}.version.txt
    """
}
