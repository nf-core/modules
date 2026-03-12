process XENIUMRANGER_RESEGMENT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/")

    output:
    tuple val(meta), path("${prefix}"), emit: outs
    tuple val("${task.process}"), val("xeniumranger"), eval("xeniumranger -V | sed -e 's/.*xenium-//'"), emit: versions_xeniumranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RESEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""

    """
    xeniumranger resegment \\
        --id="XENIUMRANGER_RESEGMENT" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    rm -rf "${prefix}"
    mv XENIUMRANGER_RESEGMENT/outs "${prefix}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/experiment.xenium"
    """
}
