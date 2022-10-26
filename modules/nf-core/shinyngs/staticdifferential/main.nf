process SHINYNGS_STATICDIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::r-shinyngs=1.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.4.1--r41hdfd78af_0':
        'quay.io/biocontainers/r-shinyngs:1.4.1--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(differential_result)                              // Differential info: contrast and differential stats
    tuple val(meta2), path(sample), path(feature_meta), path(assay_file)    // Experiment-level info

    output:
    tuple val(meta), path("*/png/volcano.png"), path("*/html/volcano.html"), path("*/html/volcano_files")   , emit: volcanos
    path "versions.yml"                                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/differential_plots.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    differential_plots.R \\
        --differential_file $differential_result \\
        --gene_metadata_file $feature_meta \\
        --reference_level $meta.reference \\
        --treatment_level $meta.target \\
        --outdir $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
