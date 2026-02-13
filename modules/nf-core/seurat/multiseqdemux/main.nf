process SEURAT_MULTISEQDEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:3.0.2--r36h0357c0b_1':
        'biocontainers/r-seurat:3.0.2--r36h0357c0b_0' }"

    input:
    tuple val(meta), path(hto_matrix)
    tuple val(meta), path(rna_matrix)
    val produce_plots
    output:
    tuple val(meta), path("*_assignment.csv")                                   ,emit: assignment
    tuple val(meta), path("*_classification.csv")                               ,emit: classification
    tuple val(meta), path("*.R_sessionInfo.log")                                ,emit: session_info
    tuple val(meta), path("*_ridge_plot.jpeg")                 ,optional:true  ,emit: ridgePlot
    tuple val(meta), path("*_feature_scatter_plot.jpeg")       ,optional:true  ,emit: featureScatter
    tuple val(meta), path("*_violin_plot.jpeg")                ,optional:true  ,emit: violingPlot
    tuple val(meta), path("*_tsne_hto_plot.jpeg")              ,optional:true  ,emit: tSNE_hto_Plot
    tuple val(meta), path("*_tsne_classification.jpeg")        ,optional:true  ,emit: tsne_classification
    path "versions.yml"                                                        ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'MULTIseq.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_assignment.csv
    touch ${prefix}_classification.csv
    touch ${prefix}_ridge_plot.jpeg
    touch ${prefix}_feature_scatter_plot.jpeg
    touch ${prefix}_violin_plot.jpeg
    touch ${prefix}_tsne_hto_plot.jpeg
    touch ${prefix}_tsne_classification.jpeg
    touch ${prefix}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
