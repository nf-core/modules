process CELLRANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    container "nf-core/cellranger:10.0.0"

    input:
    path gtf

    output:
    path "*.gtf", emit: gtf
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | sed "s/.*-//"'), emit: versions_cellranger, topic: versions

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
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${gtf.baseName}.filtered"
    """
    touch ${prefix}.gtf
    """
}
