// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bwameth=0.2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwameth:0.2.2--py_1"
    } else {
        container "quay.io/biocontainers/bwameth:0.2.2--py_1"
    }

    input:
    path fasta, stageAs: "bwameth/*"

    output:
    path "bwameth"      , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    bwameth.py index $fasta

    echo \$(bwameth.py --version 2>&1) | cut -f2 -d" " > ${software}.version.txt
    """
}
