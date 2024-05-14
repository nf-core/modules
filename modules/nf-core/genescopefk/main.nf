process GENESCOPEFK {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.2'

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GENESCOPEFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def GENESCOPE_VERSION = '380815c420f50171f9234a0fd1ff426b39829b91' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
