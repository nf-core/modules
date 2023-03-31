process CELLPOSE {
    tag "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "I did not manage to create a cellpose module in Conda that works in all OSes. Please use Docker / Singularity / Podman instead."}

    container "biocontainers/cellpose:2.1.1_cv1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*masks.tif"),   emit: mask
    path "versions.yml"                ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION='2.1.1'
    """
    cellpose \
    --image_path $image \
    --save_tif \
    --verbose \
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: $VERSION
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2.1.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_cp_masks.tif

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: $VERSION
    END_VERSIONS
    """

}
