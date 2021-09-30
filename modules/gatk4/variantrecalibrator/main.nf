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
    path resourcetbi
    val annotation
    val mode
    val rscript

    output:
    tuple val(meta), path("*.recal")   , emit: recal
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R")  , optional:true, emit: plots
    path "versions.yml"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def refCommand = ''
    def vcfCommand = ''
    def alleleSpecificCommand = ''
    def resourceCommand = ''
    def annotationCommand = ''
    def modeCommand = ''
    def rscriptCommand = ''

    def vcfList = []
    def resourceList = []
    def annotationList = []

    refCommand = fasta ? " -R ${fasta} " : ''
    vcf.each() {a -> vcfList.add(" -V " + a ) }
    vcfCommand = vcfList.join( ' ')
    if(alleleSpecific){
        alleleSpecificCommand = " --use-alle-specific-annotations"
    } else {
        alleleSpecificCommand = ''
    }
    resource.each() {a -> resourceList.add(" --resource " + a ) }
    resourceCommand = resourceList.join( ' ')
    annotation.each() {a -> annotationList.add(" -an " + a ) }
    annotationCommand = annotationList.join( ' ')
    modeCommand = mode ? " --mode ${mode} " : ''
    if(rscript){
        rscriptCommand = " --rscript-file ${prefix}.plots.R"
    } else {
        rscriptCommand = ''
    }
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
