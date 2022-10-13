process DEEPTOOLS_PLOTHEATMAP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    def container_image = "/deeptools:3.5.1--py_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: table
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotHeatmap \\
        $args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotHeatmap.pdf \\
        --outFileNameMatrix ${prefix}.plotHeatmap.mat.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotHeatmap --version | sed -e "s/plotHeatmap //g")
    END_VERSIONS
    """
}
