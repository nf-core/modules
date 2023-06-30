process DEEPTOOLS_PLOTPCA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: tab
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotPCA \\
        $args \\
        --corData $matrix \\
        --plotFile ${prefix}.plotPCA.pdf \\
        --outFileNameData ${prefix}.plotPCA.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotPCA --version | sed -e "s/plotPCA //g")
    END_VERSIONS
    """
}
