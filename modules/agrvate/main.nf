// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process AGRVATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::agrvate=1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/agrvate:1.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/agrvate:1.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}-results/${fasta.baseName}-summary.tab"), emit: summary
    path "${fasta.baseName}-results"                                                , emit: results_dir
    path "*.version.txt"                                                            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    agrvate \\
        $options.args \\
        -i $fasta

    echo \$(agrvate -v 2>&1) | sed 's/agrvate //;' > ${software}.version.txt
    """
}
