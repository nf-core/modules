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

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8a/8acde9b0bc0374d51ca50f3e49fb0f26a6fb394d277b3af5e276de8bf1786610/data'
        : 'community.wave.seqera.io/library/r-shinyngs:3.1.1--b5595b301bfbfd4d'}"

    input:
    tuple val(meta), path(sample), path(feature_meta), path(assay_files)    // Experiment-level info
    tuple val(meta2), path(contrasts), path(differential_results)           // Differential info: contrasts and differential stats
    val(contrast_stats_assay)

    output:
    tuple val(meta), path("*/data.rds"), path("*/app.R")    , emit: app
    tuple val("${task.process}"), val('shinyngs'), eval('Rscript -e "library(shinyngs); cat(as.character(packageVersion(\'shinyngs\')))"'), emit: versions_shinyngs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/make_app_from_files.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    make_app_from_files.R \\
        --sample_metadata "$sample" \\
        --feature_metadata "$feature_meta" \\
        --assay_files "${assay_files.join(',')}" \\
        --contrast_file "$contrasts" \\
        --contrast_stats_assay "$contrast_stats_assay" \\
        --differential_results "${differential_results.join(',')}" \\
        --output_dir "$prefix" \\
        $args \\
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id

    """
    mkdir -p "$prefix"
    touch "${prefix}/data.rds"
    touch "${prefix}/app.R"
    """

}
