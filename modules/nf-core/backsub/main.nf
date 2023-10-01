process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/background_subtraction:v0.3.4"

    input:
    tuple val(meta) , path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("${prefix}.backsub.ome.tif"), emit: backsub_tif
    tuple val(meta2), path("*.csv")                   , emit: markerout
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 /background_subtraction/background_sub.py \
        -o "${prefix}.backsub.ome.tif" \
        -mo markers_bs.csv \
        -r $image \
        -m $markerfile \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(python3 /background_subtraction/background_sub.py --version | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.backsub.ome.tif"
    touch "markers_bs.csv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(python3 /background_subtraction/background_sub.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
