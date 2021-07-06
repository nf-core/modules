// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ucsc-bigwigaverageoverbed=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    }

    input:
    tuple val(meta), path(bed)
    path bigwig

    output:
    tuple val(meta), path("*.tab")           , emit: tab
    path "*.version.txt"                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    # there is a bug that bigWigAverageOverBed can not handle ensembl seqlevels style.
    bigWigAverageOverBed ${options.args} $bigwig $bed ${bed.getSimpleName()}.tab

    echo \$(bigWigAverageOverBed 2>&1) | sed 's/bigWigAverageOverBed v//; s/ - Compute.*\$//' > ${software}.version.txt
    """
}
