
process GATK4_ANALYZECOVARIATES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(recal1), path(recal2), path(recal3)
    val csvout
    val ignoretimewarning

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.csv"), optional:true, emit: csv
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def inputsCommand = ''
    if(recal3) {
    inputsCommand = "-before ${recal1} -after ${recal2} -bqsr ${recal3}"
    } else if(recal2) {
    inputsCommand =  "-before ${recal1} -after ${recal2}"
    } else {
    inputsCommand = "-bqsr ${recal1}"
    }
    ignoreTimeCommand = ignoretimewarning ? "--ignore-last-modification-times" : ''
    csvCommand = csvout ? "-csv ${prefix}.csv" : ''

    """
    gatk AnalyzeCovariates  \\
        ${inputsCommand} \\
        ${ignoreTimeCommand} \\
        -plots ${prefix}.pdf \\
        ${csvCommand} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
