// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '377'

process UCSC_BED12TOBIGBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:377--h446ed27_1"
    } else {
        container "quay.io/biocontainers/ucsc-bedtobigbed:377--h446ed27_1"
    }

    input:
    tuple val(meta), path(bed)
    path  sizes

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    path "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
        ${prefix}.bigBed
    echo $VERSION > ${software}.version.txt
    """
}
