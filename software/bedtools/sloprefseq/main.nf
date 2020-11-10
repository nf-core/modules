// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [upstream: 1,
                  downstream: 10 ]

def options    = initOptions(params.options)

process BEDTOOLS_SLOPREFSEQ {
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

    input:
        tuple val(meta), path(beds)
        path  sizes

    output:
        tuple val(meta), path("*.sloprefseq.bed"), emit: bed
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def beds_files = beds.sort()
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        slopBed -i ln -s ${beds[0]} -g ln -s $sizes -l ${params.upstream} -r ${params.downstream} > ${prefix}.sloprefseq.bed
        bedtools --version | sed -e "s/Bedtools v//g" > ${software}.version.txt
        """
}
