process ILASTIK_MULTICUT {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/biocontainers/ilastik:1.4.0_cv1"

    input:
    tuple val(meta), path(h5), path(probs)
    path (ilp)

    output:
    tuple val(meta), path("*.tiff") , emit: mask
    tuple val("${task.process}"), val('ilastik'), eval("/opt/ilastik-1.4.0-Linux/run_ilastik.sh --headless --version"), emit: versions_ilastik, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "ILASTIK_MULTICUT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export MPLCONFIGDIR=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache

    /opt/ilastik-1.4.0-Linux/run_ilastik.sh \\
        --headless \\
        --readonly 1 \\
        --project=$ilp \\
        --raw_data=$h5 \\
        --probabilities=${probs} \\
        --export_source="Multicut Segmentation" \\
        --output_filename_format=${prefix}.tiff \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "ILASTIK_MULTICUT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tiff
    """
}
