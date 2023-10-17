process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::genomescope2=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py310r41hdfd78af_5':
        'biocontainers/genomescope2:2.0--py310r41hdfd78af_5' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*_linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("*_transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("*_log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("*_transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("*_model.txt")                  , emit: model
    tuple val(meta), path("*_summary.txt")                , emit: summary
    path "versions.yml"                                   , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """
}
