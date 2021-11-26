process GATK4_APPLYVQSR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
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
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
