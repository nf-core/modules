// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

        conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
        } else {
            container "quay.io/biocontainers/gatk4:4.2.0.0--0"
        }

    input:
    tuple val(meta), path(vcf) , path(tbi)

    path fasta
    path fastaidx
    path dict
    val allelespecific
    path resource

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , optional:true, emit: plots
    path "*.version.txt"               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    -R \\
    -V \\
    --use-alle-specific-annotations \\
    --resource \\
    --use-annotation \\
    --mode \\
    """
    gatk VariantRecalibrator \\
        ${refCommand} \\
        ${vcfCommand} \\
        ${alleleSpecificCommand} \\
        ${resourceCommand} \\
        ${annotationCommand} \\
        ${modeCommand} \\
        -O ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        --rscript-file ${prefix}.plots.R \\
        $options.args

    echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > ${software}.version.txt
    """
}
