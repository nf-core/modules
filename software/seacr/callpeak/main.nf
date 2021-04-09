// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '1.3'

process SEACR_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::seacr=1.3 conda-forge::r-base=4.0.2 bioconda::bedtools=2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:5bb5ed4307a8187a7f34730b00431de93688fa59-0"
    } else {
        container 'quay.io/biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:5bb5ed4307a8187a7f34730b00431de93688fa59-0'
    }

    input:
    tuple val(meta), path(bedgraph), path(ctrlbedgraph)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    SEACR_1.3.sh \\
        $bedgraph \\
        $ctrlbedgraph \\
        $options.args \\
        $prefix

    echo $VERSION > ${software}.version.txt
    """
}
