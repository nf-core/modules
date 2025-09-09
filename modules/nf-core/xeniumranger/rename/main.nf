process XENIUMRANGER_RENAME {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/xeniumranger:3.0.1"

    input:
    tuple val(meta), path(xenium_bundle)
    val(region_name)
    val(cassette_name)

    output:
    tuple val(meta), path("**/outs/**"), emit: outs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RENAME module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    xeniumranger rename \\
        --id="${prefix}" \\
        --xenium-bundle="${xenium_bundle}" \\
        --region-name="${region_name}" \\
        --cassette-name="${cassette_name}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "XENIUMRANGER_RENAME module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    touch "${prefix}/outs/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xeniumranger: \$(xeniumranger -V | sed -e "s/xeniumranger-/- /g")
    END_VERSIONS
    """
}
