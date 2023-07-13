process HICEXPLORER_HICPLOTMATRIX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 numpy=1.23.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1':
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix), path(additional_files)

    output:
    tuple val(meta), path("${prefix}.${postfix}")   , emit:plots
    path("versions.yml")                            , emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    postfix = task.ext.postfix ?: 'png'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    hicPlotMatrix \\
        -m $matrix \\
        $args \\
        -o ${prefix}.${postfix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotMatrix --version 2>&1 | sed 's/hicPlotMatrix //')
    END_VERSIONS
    """
}
