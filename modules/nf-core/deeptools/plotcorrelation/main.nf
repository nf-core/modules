process DEEPTOOLS_PLOTCORRELATION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(matrix)
    val(method)
    val(plot_type)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: matrix
    tuple val("${task.process}"), val('deeptools'), eval('plotCorrelation --version | sed "s/plotCorrelation //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def resolved_method = method ?: 'spearman'
    def resolved_plot_type = plot_type ?: 'heatmap'
    """
    plotCorrelation \\
        $args \\
        --corData $matrix \\
        --corMethod $resolved_method \\
        --whatToPlot $resolved_plot_type \\
        --plotFile ${prefix}.plotCorrelation.pdf \\
        --outFileCorMatrix ${prefix}.plotCorrelation.mat.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotCorrelation.pdf
    touch ${prefix}.plotCorrelation.mat.tab
    """
}
