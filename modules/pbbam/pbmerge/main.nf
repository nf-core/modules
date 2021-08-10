// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PBBAM_PBMERGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pbbam=1.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pbbam:1.6.0--h058f120_1"
    } else {
        container "quay.io/biocontainers/pbbam:1.6.0--h058f120_1"
    }

    input:
    tuple val(meta), path("*.bam")

    output:
    tuple val(meta), path("${meta.id}.ccs.bam"), path("*.bam.pbi"), emit: bam
    path "*.version.txt"                             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def out_bam  = "${prefix}.ccs.bam"
    """
    pbmerge -o ${out_bam} $options.args *.bam
    pbindex ${out_bam}
    echo \$(pbmerge --version 2>&1) > ${software}.version.txt
    """
}
