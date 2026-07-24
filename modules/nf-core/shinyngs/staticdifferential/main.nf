process SHINYNGS_STATICDIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5f79ef0afe3e3831496c61c81aeda56312f49ac324dc378d3312af7acae2ec6/data'
        : 'community.wave.seqera.io/library/r-shinyngs:3.2.1--d43071e62bc500d3'}"

    input:
    tuple val(meta), path(differential_result)                              // Differential info: contrast and differential stats
    tuple val(meta2), path(sample), path(feature_meta), path(assay_file)    // Experiment-level info

    output:
    tuple val(meta), path("*/png/volcano.png")      , emit: volcanos_png
    tuple val(meta), path("*/html/volcano.html")    , emit: volcanos_html, optional: true
    tuple val("${task.process}"), val('shinyngs'), eval('Rscript -e "library(shinyngs); cat(as.character(packageVersion(\'shinyngs\')))"'), emit: versions_shinyngs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // For full list of available args see
    // https://github.com/pinin4fjords/shinyngs/blob/develop/exec/differential_plots.R
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    differential_plots.R \\
        --differential_file "$differential_result" \\
        --feature_metadata "$feature_meta" \\
        --outdir "$prefix" \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p $prefix/png && mkdir $prefix/html
    touch $prefix/png/volcano.png
    touch $prefix/html/volcano.html
    """
}
