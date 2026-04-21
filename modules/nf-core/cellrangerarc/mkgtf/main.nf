process CELLRANGERARC_MKGTF {
    tag "$meta.id"
    label 'process_low'

    container "nf-core/cellranger-arc:2.0.2"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    tuple val("${task.process}"), val('cellrangerarc'), eval("cellranger-arc --version 2>&1 | sed 's/cellranger-arc cellranger-arc-//'"), emit: versions_cellrangerarc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.filtered"
    """
    cellranger-arc \\
        mkgtf \\
        ${gtf} \\
        ${prefix}.gtf \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    """
    touch ${prefix}.gtf
    """
}
