process BACKSUB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/background_subtraction:v0.5.1"

    input:
    tuple val(meta) , path(image)
    tuple val(meta2), path(markerfile)

    output:
    tuple val(meta), path("*.ome.tif"), emit: backsub_tif
    tuple val(meta2), path("*.csv")   , emit: markerout
    tuple val("${task.process}"), val('backsub'), eval('backsub --version'), emit: versions_backsub, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$image" == "${prefix}.ome.tif") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    backsub \
        --input $image \
        --markers $markerfile \
        --output "${prefix}.ome.tif" \
        --marker-output "${prefix}.csv" \
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ome.tif"
    touch "${prefix}.csv"
    """
}
