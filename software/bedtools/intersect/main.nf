// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

def options    = initOptions(params.options)

process BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda     (params.enable_conda ? "bioconda::bedtools =2.29.2" : null)
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"

    input:
        tuple val(meta), path(beds)

    output:
        tuple val(meta), path("*.intersect.bed"), emit: intersect
        path  "*.version.txt", emit: version

    script: // TODO change script to account for multiple possible intersections
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        bedtools intersect -a ${beds[0]} -b ${beds[1]} ${options.args} > ${prefix}.intersect.bed
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
}
