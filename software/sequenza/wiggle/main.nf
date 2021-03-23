// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQUENZA_WIGGLE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py38h6ed170a_2"
    } else {
        container "quay.io/biocontainers/sequenza-utils:3.0.0--py38h6ed170a_2"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig.gz"), emit: wig
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    sequenza-utils \\
        gc_wiggle \\
        $options.args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz

    echo \$(sequenza-utils --version 2>&1) | sed 's/^.*sequenza-utils //; s/Using.*\$//' > ${software}.version.txt
    """
}
