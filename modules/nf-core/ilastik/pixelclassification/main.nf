process ILASTIK_PIXELCLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "ILASTIK_PIXELCLASSIFICATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "labsyspharm/mcmicro-ilastik:1.6.1"

    input:
    tuple val(meta), path(h5)
    tuple val(meta2), path(ilp)

    output:
    tuple val(meta), path("*.${suffix}") , emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "h5"

    """
    /ilastik-release/run_ilastik.sh \\
        --headless \\
        --readonly 1 \\
        --project=$ilp \\
        --output_filename_format=${prefix}.${suffix}
        $args \\
        $h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(/ilastik-release/run_ilastik.sh --headless --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "h5"

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik:: \$(/ilastik-release/run_ilastik.sh --headless --version)
    END_VERSIONS
    """
}
