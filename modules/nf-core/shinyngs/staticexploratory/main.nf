process SHINYNGS_STATICEXPLORATORY {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::r-shinyngs=1.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs%3A1.3.2--r41hdfd78af_1':
        'quay.io/biocontainers/r-shinyngs:1.3.2--r41hdfd78af_1' }"

    input:
    tuple val(meta), path(sample), path(feature_meta), path(assay_files) 

    output:
    tuple val(meta), path("*/png/boxplot.png"), path("*/html/boxplot.html"), path("*/html/boxplot_files")                         , emit: boxplots
    tuple val(meta), path("*/png/density.png"), path("*/html/density.html"), path("*/html/density_files")                         , emit: densities
    tuple val(meta), path("*/png/pca2d.png"), path("*/html/pca2d.html"), path("*/html/pca2d_files")                               , emit: pca2d
    tuple val(meta), path("*/png/pca3d.png"), path("*/html/pca3d.html"), path("*/html/pca3d_files")                               , optional: true, emit: pca3d
    tuple val(meta), path("*/png/mad_correlation.png"), path("*/html/mad_correlation.html"), path("*/html/mad_correlation_files") , emit: mad
    tuple val(meta), path("*/png/sample_dendrogram.png")                                                                          , emit: dendro
    path "versions.yml"                                                                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/exploratory_plots.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    exploratory_plots.R \\
        --sample_metadata $sample \\
        --feature_metadata $feature_meta \\
        --assay_files ${assay_files.join(',')} \\
        --contrast_variable ${meta.id} \\
        --outdir $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
