process HICEXPLORER_HICCORRECTMATRIX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2 numpy=1.23.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1':
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${prefix}.${postfix}")   , emit:corrected
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    postfix = task.ext.postfix ?: 'h5'
    prefix = task.ext.prefix ?: "${meta.id}_corrected"
    """
    hicCorrectMatrix correct \\
        -m $matrix \\
        $args \\
        -o ${prefix}.${postfix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicCorrectMatrix --version 2>&1 | sed 's/hicCorrectMatrix //')
    END_VERSIONS
    """
}
