process GENESCOPEFK {
    tag "$meta.id"
    label 'process_low'

    if (params.enable_conda) {
        error "Conda environments cannot be used when using the GeneScope tool. Please use docker or singularity containers."
    }
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.0'

    input:
    tuple val(meta), path(fastk_histex_histogram)

    output:
    tuple val(meta), path("*_linear_plot.png")            , emit: linear_plot
    tuple val(meta), path("*_log_plot.png")               , emit: log_plot
    tuple val(meta), path("*_model.txt")                  , emit: model
    tuple val(meta), path("*_summary.txt")                , emit: summary
    tuple val(meta), path("*_transformed_linear_plot.png"), emit: transformed_linear_plot
    tuple val(meta), path("*_transformed_log_plot.png")   , emit: transformed_log_plot
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def GENESCOPE_VERSION = '380815c420f50171f9234a0fd1ff426b39829b91' // Git commit id is used instead of GeneScopeFK.R -v as software is not release versioned.
    """
    GeneScopeFK.R \\
        $args \\
        --input $fastk_histex_histogram \\
        --output . \\
        --name_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genescope: $GENESCOPE_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
