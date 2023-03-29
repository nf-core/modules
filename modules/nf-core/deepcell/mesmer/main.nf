process MESMER {
    tag "$meta.id"
    label 'process_medium'

    // Setting up the container to latest all the time because why not.
    container "vanvalenlab/deepcell-applications:0.4.0"

    // Mesmer, requieres one image to segment and the mpp(microns per pixel)
    input:
    tuple val(meta), path(img)

    // Output a .tif image, don't touch versions
    output:
    tuple val(meta), path("mask.tif"), emit: mask
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.4.0"

    """
    python /usr/src/app/run_app.py mesmer \
        --squeeze \
        --nuclear-image $img \
        --membrane-image $img \
        --output-directory . \
        --output-name mask.tif \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepcell_mesmer:: $VERSION
    END_VERSIONS
    """
}
