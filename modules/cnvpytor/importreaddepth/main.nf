process CNVPYTOR_IMPORTREADDEPTH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::cnvpytor=1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:A1.0--py39h6a678da_2':
        'quay.io/biocontainers/cnvpytor:1.0--py39h6a678da_2' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.pytor")	, emit: pytor
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cnvpytor \\
        -root ${prefix}.pytor \\
        -rd $bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/^.*cnvpytor //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
