// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAPHOPIASCCMEC {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::staphopia-sccmec=1.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/staphopia-sccmec:1.0.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    staphopia-sccmec --assembly $fasta $options.args > ${prefix}.tsv

    echo \$(staphopia-sccmec --version 2>&1) | sed 's/^.*staphopia-sccmec //' > ${software}.version.txt
    """
}
