// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

        conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
        } else {
            container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
        }

    input:
    tuple val(meta), path(vcf) , path(tbi)
    path fasta
    path fai
    path dict
    val allelespecific
    tuple path(resvcfs), path(restbis), val(reslabels)
    val annotation
    val mode
    val create_rscript

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , optional:true, emit: plots
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    refCommand = fasta ? "-R ${fasta} " : ''
    vcfCommand = '-V ' + vcf.join( ' -V ')
    alleleSpecificCommand = allelespecific ? '-AS' : ''
    resourceCommand = '--resource:' + reslabels.join( ' --resource:')
    annotationCommand = '-an ' + annotation.join( ' -an ')
    modeCommand = mode ? "--mode ${mode} " : 'SNP'
    rscriptCommand = create_rscript ? "--rscript-file ${prefix}.plots.R" : ''

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
        ${rscriptCommand}\\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
