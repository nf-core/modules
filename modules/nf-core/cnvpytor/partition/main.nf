process CNVPYTOR_PARTITION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.3.1--pyhdfd78af_1':
        'biocontainers/cnvpytor:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(pytor)
    val bin_sizes

    output:
    tuple val(meta), path("${prefix}.pytor"), emit: pytor
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bins = bin_sizes       ?: '1000'
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(pwd)/.mplconfig

    cnvpytor \\
        -root ${pytor} \\
        -partition ${bins} \\
        ${args}

    cp ${pytor} ${prefix}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version | sed -n 's/.*CNVpytor \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
