process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe5ecaf5a34b7345080a9b54be83711b9a9732fdcbfd1809338da6b277bb9dca/data':
        'community.wave.seqera.io/library/genomescope2:2.1.0--a4c756d0a4552c53' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("${prefix}_linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("${prefix}_transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("${prefix}_log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("${prefix}_transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("${prefix}_model.txt")                  , emit: model
    tuple val(meta), path("${prefix}_summary.txt")                , emit: summary
    tuple val(meta), path("${prefix}_lookup_table.txt")           , emit: lookup_table, optional: true
    tuple val(meta), path("${prefix}_fitted_hist.png")            , emit: fitted_histogram_png, optional: true
    tuple val(meta), path("*.json")                               , emit: json_report, optional: true
    tuple val("${task.process}"), val('genomescope2'), eval('genomescope2 -v | sed "s/GenomeScope //"'), emit: versions_genomescope2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomescope2 \\
        --input $histogram \\
        $args \\
        --output . \\
        --name_prefix $prefix

    test -f "fitted_hist.png" && mv fitted_hist.png ${prefix}_fitted_hist.png
    test -f "lookup_table.txt" && mv lookup_table.txt ${prefix}_lookup_table.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_linear_plot.png
    touch ${prefix}_transformed_linear_plot.png
    touch ${prefix}_log_plot.png
    touch ${prefix}_transformed_log_plot.png
    touch ${prefix}_model.txt
    touch ${prefix}_summary.txt
    touch ${prefix}_fitted_hist.png
    touch ${prefix}_lookup_table.txt

    echo ${args}
    """
}
