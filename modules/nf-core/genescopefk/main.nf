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
    tuple val(meta), path("*_genescopefk.log")            , emit: log
    tuple val(meta), env('KMERCOV')                       , emit: kmer_cov
    tuple val("${task.process}"), val('genescopefk'), eval("R --version | sed '1!d; s/.*version //; s/ .*//'"), emit: versions_genescopefk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GENESCOPEFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    GeneScopeFK.R \\
        $args \\
        --input $fastk_histex_histogram \\
        --output . \\
        --name_prefix ${prefix} > ${prefix}_genescopefk.log

    printf -v KMERCOV "%.2f" \$( grep "^kmercov" *_model.txt | cut -d" " -f2 )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_linear_plot.png"
    touch "${prefix}_log_plot.png"
    touch "${prefix}_model.txt"
    touch "${prefix}_summary.txt"
    touch "${prefix}_transformed_linear_plot.png"
    touch "${prefix}_transformed_log_plot.png"
    touch "${prefix}_genescopefk.log"

    printf -v KMERCOV "%.2f" \$( grep "^kmercov" *_model.txt | cut -d" " -f2 )
    """

}
