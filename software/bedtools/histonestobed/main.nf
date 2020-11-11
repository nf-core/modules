// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BEDTOOLS_HISTONESTOBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bedtools =2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0 "
    } else {
        container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    }
    // In the case of use for cell lines, the second input can also be replaced with the tuple (val) format.
    input:
        tuple val(meta), path(bams)

    output:
        tuple val(meta), path("*.bed"), emit: bed
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        bamToBed -i ${bams[0]} | sortBed -i > ${prefix}.bed 
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
}