process CELLRANGER_MKGTF {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/nf-core/cellranger:10.0.0"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | sed "s/.*-//"'), emit: versions_cellranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKGTF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.filtered"
    if ("${gtf}" == "${prefix}.gtf") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
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
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    if ("${gtf}" == "${prefix}.gtf") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.gtf
    """
}
