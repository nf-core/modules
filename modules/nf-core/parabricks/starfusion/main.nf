process PARABRICKS_STARFUSION {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta), path(chimeric_junction)
    tuple val(meta1), path(genome_lib_dir)

    output:
    tuple val(meta), path("fusion_predictions.tsv"),                emit: fusions
    tuple val(meta), path("fusion_predictions.abridged.tsv"),       emit: abridged
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        starfusion \\
        --chimeric-junction ${chimeric_junction} \\
        --genome-lib-dir ${genome_lib_dir} \\
        --output-dir ${prefix}_starfusion \\
        ${num_gpus} \\
        ${args}

    mv ${prefix}_starfusion/fusion_predictions.tsv .
    mv ${prefix}_starfusion/fusion_predictions.abridged.tsv .
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    touch fusion_predictions.tsv
    touch fusion_predictions.abridged.tsv
    """
}
