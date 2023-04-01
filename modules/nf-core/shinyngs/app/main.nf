
process SHINYNGS_APP {
    tag "$meta.id"
    label 'process_single'

    // To be able to pass the necessary secrets for shinyapps.io deployment,
    // this process must be configured by placing a statement like the
    // following in the nextflow.config:
    //
    // withName: SHINYNGS_APP {
    //     secret = [ 'SHINYAPPS_TOKEN', 'SHINYAPPS_SECRET' ]
    // }
    //
    // Those values must then be set in your Nextflow secrets.

    conda "bioconda::r-shinyngs=1.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.6.0--r42hdfd78af_1':
        'quay.io/biocontainers/r-shinyngs:1.6.0--r42hdfd78af_1' }"

    input:
    tuple val(meta), path(sample), path(feature_meta), path(assay_files)    // Experiment-level info
    tuple val(meta2), path(contrasts), path(differential_results)           // Differential info: contrasts and differential stats
    val(contrast_stats_assay)

    output:
    tuple val(meta), path("*/data.rds"), path("*/app.R")    , emit: app
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/make_app_from_files.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    make_app_from_files.R \\
        --sample_metadata $sample \\
        --feature_metadata $feature_meta \\
        --assay_files ${assay_files.join(',')} \\
        --contrast_file $contrasts \\
        --contrast_stats_assay $contrast_stats_assay \\
        --differential_results ${differential_results.join(',')} \\
        --output_dir $prefix \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
