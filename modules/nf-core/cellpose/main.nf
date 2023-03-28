process CELLPOSE {
    tag "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "I did not manage to create a cellpose module in Conda that works in all OSes. Please use Docker / Singularity / Podman instead."
    }

    container "biocontainers/cellpose:2.1.1_cv1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*mask.tif"), emit: mcquant

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cellpose \
    --image_path $image \
    --channel_axis=0 \
    --no_npy \
    --save_tif \
    --verbose
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2.1.1_cv1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
}
