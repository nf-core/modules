process DEEPTOOLS_PLOTPCA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('deeptools'), eval('plotPCA --version | sed "s/plotPCA //g"') , emit: versions_deeptools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotPCA.pdf
    touch ${prefix}.plotPCA.tab
    """
}
