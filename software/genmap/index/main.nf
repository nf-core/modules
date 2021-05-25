// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENMAP_INDEX {
    tag '$fasta'
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::genmap=1.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/genmap:1.3.0--h1b792b2_1"
    } else {
        container "quay.io/biocontainers/genmap:1.3.0--h1b792b2_1"
    }

    input:
    path fasta

    output:
    path "genmap"       , emit: index
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    genmap \\
        index \\
        -F $fasta \\
        -I genmap

    echo \$(genmap --version 2>&1) | sed 's/GenMap version: //; s/SeqAn.*\$//' > ${software}.version.txt
    """
}
