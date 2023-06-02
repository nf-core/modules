process SPACERANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    // TODO push to nf-core docker
    container "ghcr.io/grst/spaceranger:2.1.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "SPACERANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path gtf

    output:
    path "*.gtf"         , emit: gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}.filtered"
    """
    spaceranger \\
        mkgtf \\
        $gtf \\
        ${prefix}.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """
}
