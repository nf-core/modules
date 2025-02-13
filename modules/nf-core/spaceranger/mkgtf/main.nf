process SPACERANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    container "nf-core/spaceranger:3.1.3"

    input:
    path gtf

    output:
    path "*.gtf"         , emit: gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
