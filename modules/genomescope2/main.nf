process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::genomescope2=2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py310r41hdfd78af_5':
        'quay.io/biocontainers/genomescope2:2.0--py310r41hdfd78af_5' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("$prefix-linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("$prefix-transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("$prefix-log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("$prefix-transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("$prefix-model.txt")                  , emit: model
    tuple val(meta), path("$prefix-summary.txt")                , emit: summary
    path "versions.yml"                                         , emit: versions

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

    ls -l

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """
}
