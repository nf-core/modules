process XENIUMRANGER_RENAME {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/"), val(region_name), val(cassette_name)

    output:
    tuple val(meta), path("${prefix}"), emit: outs
    tuple val("${task.process}"), val("xeniumranger"), eval("xeniumranger -V | sed -e 's/.*xenium-//'"), emit: versions_xeniumranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RENAME module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    rm -rf "${prefix}"

    xeniumranger rename \\
        --id="XENIUMRANGER_RENAME" \\
        --xenium-bundle="${xenium_bundle}" \\
        --region-name="${region_name}" \\
        --cassette-name="${cassette_name}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    mv XENIUMRANGER_RENAME/outs "${prefix}"
    """

    stub:

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "${prefix}"
    touch "${prefix}/experiment.xenium"
    """
}
