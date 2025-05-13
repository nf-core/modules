process MCSTAGING_PHENOIMAGER2MC {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/phenoimager2mc:v0.2.2"

    input:
    tuple val(meta) , path(tiles, stageAs: "tiles/*")

    output:
    tuple val(meta), path("*.tif"), emit: tif
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Phenoimager2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /phenoimager2mc/scripts/phenoimager2mc.py \
        -i ${tiles} \
        -o "${prefix}.tif" \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phenoimager2mc: \$(python /phenoimager2mc/scripts/phenoimager2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Phenoimager2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch input
    touch "${prefix}.tif"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phenoimager2mc: \$(python /phenoimager2mc/scripts/phenoimager2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
