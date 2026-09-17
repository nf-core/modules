process METAPROVIZ_INTERNALSTANDARD {
    tag "$meta.id"
    label 'process_single'

    // needs to be changed to "quay.io/nf-core/metaproviz:0.0.1". We need need a core member to do this
    container "ghcr.io/saezlab/metaproviz:0.0.1"

    input:
    tuple val(meta), path(se_rds), path(data_matrix), path(feature_matrix), path(sample_matrix)

    output:
    tuple val(meta), path("*.cv.tsv")            , emit: cv_table
    tuple val(meta), path("*.high_var.txt")      , emit: high_var
    tuple val(meta), path("*.condition_cv.tsv")  , emit: condition_cv
    tuple val(meta), path("*.plots.rds")         , emit: plots
    tuple val(meta), path("*.report.html")       , emit: report
    tuple val(meta), path("*.log")               , emit: log
    path "versions.yml"                                            , emit: versions_internalstandard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "METAPROVIZ_INTERNALSTANDARD module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    template 'internal_standard.R'

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "METAPROVIZ_INTERNALSTANDARD module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.cv.tsv ${prefix}.high_var.txt \\
          ${prefix}.condition_cv.tsv ${prefix}.plots.rds \\
          ${prefix}.report.html ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.1
        metaproviz: $VERSION
    END_VERSIONS
    """
}
