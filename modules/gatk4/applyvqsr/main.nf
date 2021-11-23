// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_APPLYVQSR {
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
    tuple val(meta), path(vcf), path(tbi), path(recal), path(recalidx), path(tranches)
    path fasta
    path fai
    path dict
    val allelespecific
    val truthsensitivity
    val mode

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    path "versions.yml"                   , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    refCommand = fasta ? "-R ${fasta} " : ''
    alleleSpecificCommand = allelespecific ? '-AS' : ''
    truthSensitivityCommand = truthsensitivity ? "--truth-sensitivity-filter-level ${truthsensitivity}" : ''
    modeCommand = mode ? "--mode ${mode} " : 'SNP'
    """
    gatk ApplyVQSR \\
        ${refCommand} \\
        -V ${vcf} \\
        -O ${prefix}.vcf.gz \\
        ${alleleSpecificCommand} \\
        ${truthSensitivityCommand} \\
        --tranches-file $tranches \\
        --recal-file $recal \\
        ${modeCommand} \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
