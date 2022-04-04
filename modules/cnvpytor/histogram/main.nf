process CNVPYTOR_HISTOGRAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cnvpytor=1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.0--py39h6a678da_2':
        'quay.io/biocontainers/cnvpytor:1.0--py39h6a678da_2' }"

    input:
    tuple val(meta), path(pytor)

    output:
    tuple val(meta), path("${pytor.baseName}.pytor")	, emit: pytor
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '1000'
    """
    cnvpytor \\
        -root $pytor \\
        -his $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/^.*pyCNVnator //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    """
    touch ${pytor.baseName}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/^.*pyCNVnator //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
