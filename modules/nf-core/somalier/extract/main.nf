
process SOMALIER_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/somalier:0.2.19--h0c29559_0':
        'biocontainers/somalier:0.2.19--h0c29559_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(sites)

    output:
    tuple val(meta), path("*.somalier") , emit: extract
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    somalier extract \\
        --sites ${sites} \\
        -f ${fasta} \\
        ${input} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.somalier

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
    """
}
