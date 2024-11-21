process SKA_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ska:1.0--he860b03_0':
        'biocontainers/ska:1.0--he860b03_0' }"

    input:
    tuple val(meta), path(fastas), path(fasta_list)

    output:
    tuple val(meta), path("*.skf"), emit: skf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def arg_list = fasta_list ? "-f ${fasta_list}" : ''
    """
    ska fasta \\
        ${args} \\
        ${fastas} \\
        ${arg_list} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska version | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.skf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska version | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
