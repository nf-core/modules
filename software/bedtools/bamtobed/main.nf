// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BEDTOOLS_BAMTOBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda     (params.enable_conda ? "bioconda::bedtools=2.29.2" : null)
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"

    input:
    tuple val(meta), path(sizes), path(bam), path(bai)
   
    output:
    tuple val(meta), path(sizes), path("*.bed12"), emit: bed12
    path "*.version.txt"                         , emit: version

    script:
    """
    bedtools \\
        bamtobed \\
        -bed12 \\
        -cigar \\
        -i ${bam[0]} \\
        | bedtools sort > ${meta.id}.bed12
    bedtools --version | sed -e "s/bedtools v//g" > bedtools.version.txt
    """
}
