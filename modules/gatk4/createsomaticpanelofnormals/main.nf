// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_CREATESOMATICPANELOFNORMALS {
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
    tuple val(meta), path(vcf)

    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputsList = []
    def inputsCommand = ''
    vcf.each() {a -> inputsList.add(" -V " + a ) }
    inputsCommand = inputsList.join(' ')

    """
    gatk CreateSomaticPanelOfNormals\\
        -R $fasta
        ${inputsCommand} \\
        -O ${prefix}.vcf.gz \\
        $options.args

    echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > ${software}.version.txt
    """
}
