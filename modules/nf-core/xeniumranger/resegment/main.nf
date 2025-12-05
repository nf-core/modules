process XENIUMRANGER_RESEGMENT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/")

    output:
    tuple val(meta), path("${prefix}"), emit: outs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RESEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = (task.ext.args ?: "").trim()
    def args_block = args ? "\\\n        ${args}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    rm -rf "${prefix}"

    xeniumranger resegment \\
        --id="XENIUMRANGER_RESEGMENT" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}${args_block}

    if [ -d "XENIUMRANGER_RESEGMENT/outs" ]; then
        rm -rf "${prefix}"
        mv XENIUMRANGER_RESEGMENT/outs "${prefix}"
    else
        mkdir -p "${prefix}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RESEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = (task.ext.args ?: "").trim()
    def args_block = args ? "\\\n        ${args}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rm -rf "${prefix}"

    xeniumranger resegment \\
        --id="XENIUMRANGER_RESEGMENT" \\
        --xenium-bundle="${xenium_bundle}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}${args_block} \\
        --dry

    if [ -d "XENIUMRANGER_RESEGMENT/outs" ]; then
        rm -rf "${prefix}"
        mv XENIUMRANGER_RESEGMENT/outs "${prefix}"
    else
        mkdir -p "${prefix}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """
}
