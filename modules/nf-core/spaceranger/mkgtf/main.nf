process SPACERANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    container "nf-core/spaceranger:9c5e7dc93c32448e"

    input:
    path(gtf)

    output:
    path("*.gtf"), emit: gtf
    tuple val("${task.process}"), val('spaceranger'), eval('spaceranger -V | sed  "s/spaceranger spaceranger-//"'), emit: versions_spaceranger, topic: versions

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
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${gtf.baseName}.filtered"
    """
    touch ${prefix}.gtf
    """
}
