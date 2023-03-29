process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/background_subtraction:v0.3.4"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("*.ome.tif"), emit: backsub_tif
    tuple val(meta2), path("*.csv")  , emit: markerout
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 'echo python3 /background_subtraction/background_sub.py --version'
    """
    python3 /background_subtraction/background_sub.py \
        -o ./img_backsub.ome.tif \
        -mo ./markerout.csv \
        -r $image \
        -m $markerfile \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: ${VERSION:1}
    END_VERSIONS
    """
}
