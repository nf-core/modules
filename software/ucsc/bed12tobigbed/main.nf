// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process UCSC_BED12TOBIGBED {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'ucsc', publish_id:'') }

    conda     (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    container "quay.io/biocontainers/ucsc-bedtobigbed:377--h446ed27_1"
    if ({ task.exitStatus in [255] }) { errorStrategy = 'ignore' }

    input:
    tuple val(sample), path(sizes), path(bed12)

    output:
    tuple val(sample), path(sizes), path("*.bigBed"), emit: bigbed

    script:
    """
    bedToBigBed \\
        $bed12 \\
        $sizes \\
        ${sample}.bigBed
    """
}
