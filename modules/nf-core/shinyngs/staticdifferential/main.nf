process SHINYNGS_STATICDIFFERENTIAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b43cb5832431a3d581f3abdca0f5d45fa814b8a4b2c4e49f12296d2e07629fd/data' :
        'community.wave.seqera.io/library/r-shinyngs:2.0.0--feeacac14bcf4bb0' }"

    input:
    tuple val(meta), path(differential_result)                              // Differential info: contrast and differential stats
    tuple val(meta2), path(sample), path(feature_meta), path(assay_file)    // Experiment-level info

    output:
    tuple val(meta), path("*/png/volcano.png")      , emit: volcanos_png
    tuple val(meta), path("*/html/volcano.html")    , emit: volcanos_html, optional: true
    path "versions.yml"                             , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p $prefix/png && mkdir $prefix/html
    touch $prefix/png/volcano.png
    touch $prefix/html/volcano.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-shinyngs: \$(Rscript -e "library(shinyngs); cat(as.character(packageVersion('shinyngs')))")
    END_VERSIONS
    """
}
