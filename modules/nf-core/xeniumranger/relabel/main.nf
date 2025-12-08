process XENIUMRANGER_RELABEL {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:4.0"

    input:
    tuple val(meta), path(xenium_bundle), path(gene_panel)
    output:
    tuple val(meta), path("${prefix}"), emit: outs
    tuple val("${task.process}"), val("xeniumranger"), eval("xeniumranger -V | sed -e 's/xeniumranger-/- /g'"), emit: versions_xeniumranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RELABEL module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = (task.ext.args ?: '').trim()
    def args_block = args ? "\\\n        ${args}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rm -rf "${prefix}"

    xeniumranger relabel \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --panel="${gene_panel}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}${args_block}

    if [ ! -d "${prefix}" ]; then
        mkdir -p "${prefix}"
    fi
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RELABEL module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = (task.ext.args ?: '').trim()
    def args_block = args ? "\\\n        ${args}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rm -rf "${prefix}"

    xeniumranger relabel \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --panel="${gene_panel}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}${args_block} \\
        --dry

    if [ ! -d "${prefix}" ]; then
        mkdir -p "${prefix}"
    fi
    """
}
