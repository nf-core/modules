process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/background_subtraction:v0.4.1"

    input:
    tuple val(meta) , path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("*.ome.tif"), emit: backsub_tif
    tuple val(meta2), path("*.csv")   , emit: markerout
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Backsub module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$image" == "${prefix}.ome.tif") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    python3 /background_subtraction/background_sub.py \
        -o "${prefix}.ome.tif" \
        -mo "${prefix}.csv" \
        -r $image \
        -m $markerfile \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(python3 /background_subtraction/background_sub.py --version | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Backsub module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ome.tif"
    touch "${prefix}.csv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(python3 /background_subtraction/background_sub.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
