process SHINYNGS_STATICEXPLORATORY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-shinyngs=1.5.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.5.5--r42hdfd78af_0':
        'quay.io/biocontainers/r-shinyngs:1.5.5--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(sample), path(feature_meta), path(assay_files)

    output:
    tuple val(meta), path("*/png/boxplot.png")                  , emit: boxplots_png
    tuple val(meta), path("*/html/boxplot.html")                , emit: boxplots_html, optional: true
    tuple val(meta), path("*/png/density.png")                  , emit: densities_png
    tuple val(meta), path("*/html/density.html")                , emit: densities_html, optional: true
    tuple val(meta), path("*/png/pca2d.png")                    , emit: pca2d_png
    tuple val(meta), path("*/html/pca2d.html")                  , emit: pca2d_html, optional: true
    tuple val(meta), path("*/png/pca3d.png")                    , emit: pca3d_png
    tuple val(meta), path("*/html/pca3d.html")                  , emit: pca3d_html, optional: true
    tuple val(meta), path("*/png/mad_correlation.png")          , emit: mad_png
    tuple val(meta), path("*/html/mad_correlation.html")        , emit: mad_html, optional: true
    tuple val(meta), path("*/png/sample_dendrogram.png")        , emit: dendro
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/exploratory_plots.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    exploratory_plots.R \\
        --sample_metadata "$sample" \\
        --feature_metadata "$feature_meta" \\
        --assay_files "${assay_files.join(',')}" \\
        --contrast_variable "${meta.id}" \\
        --outdir "$prefix" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
