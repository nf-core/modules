process CNVPYTOR_HISTOGRAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.2.1--pyhdfd78af_0':
        'biocontainers/cnvpytor:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pytor)
    val bin_sizes


    output:
    tuple val(meta), path("${pytor.baseName}.pytor")	, emit: pytor
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bins = bin_sizes ?: '1000'
    """
    cnvpytor \\
        -root $pytor \\
        -his $bins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\)\$/\\1/p' )
    END_VERSIONS
    """

    stub:
    """
    touch ${pytor.baseName}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\)\$/\\1/p' )
    END_VERSIONS
    """
}
