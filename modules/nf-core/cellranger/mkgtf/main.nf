process CELLRANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    container "docker.io/nfcore/cellranger:7.1.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path gtf

    output:
    path "*.filtered.gtf", emit: gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellranger \\
        mkgtf \\
        $gtf \\
        ${gtf.baseName}.filtered.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
