process GAWK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(input)
    path(program_file)

    output:
    tuple val(meta), path("${prefix}")  , emit: output
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: '' // args is used for the main arguments of the tool
    def args2 = task.ext.args2 ?: '' // args2 is used to specify a program when no program file has been given
    prefix = task.ext.prefix ?: "${meta.id}.${input.getExtension}" // This default is up to debate, please post about it on Slack or open a PR with your suggested change

    program = program_file ? "-f ${program_file}" : "${args2}"

    """
    awk \\
        ${args} \\
        ${program} \\
        ${input} \\
        > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
