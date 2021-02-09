// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_FILTERVARIANTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/1.3.1--h089eab3_0"
    } else {
        container "quay.io/biocontainers/1.3.1--h089eab3_0"
    }

    input:
    tuple val(meta), path(variant_tsv)

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: filtered_tsv
    path "*.version.txt"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // is the input a file of filenames (one per line)?
    def fofn     = options.fofn ? "-f " : ""
    """
    ivar filtervariants \\
        $options.args \\
        -p $prefix \\
        ${fofn}${variant_tsv}

    ivar version | head -n1 2>&1 | sed 's/^.*iVar version //g' > ${software}.version.txt
    """
}
