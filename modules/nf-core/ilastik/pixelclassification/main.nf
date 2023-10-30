process ILASTIK_PIXELCLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/biocontainers/ilastik:1.4.0_cv1"

    input:
    tuple val(meta), path(input_img)
    tuple val(meta2), path(ilp)

    output:
    tuple val(meta), path("*.${suffix}") , emit: output
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "ILASTIK_PIXELCLASSIFICATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "h5"

    """
    /opt/ilastik-1.4.0-Linux/run_ilastik.sh \\
        --headless \\
        --readonly 1 \\
        --project=$ilp \\
        --output_filename_format=${prefix}.${suffix} \\
        $args \\
        $input_img

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(/opt/ilastik-1.4.0-Linux/run_ilastik.sh --headless --version)
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "ILASTIK_PIXELCLASSIFICATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "h5"

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik:: \$(/opt/ilastik-1.4.0-Linux/run_ilastik.sh --headless --version)
    END_VERSIONS
    """
}
