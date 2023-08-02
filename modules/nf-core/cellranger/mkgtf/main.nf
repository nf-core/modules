process CELLRANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    container "nf-core/cellranger:7.1.0"

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
        error "CELLRANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gtf.baseName}.filtered"
    """
    cellranger \\
        mkgtf \\
        $gtf \\
        ${prefix}.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
