process COREOGRAPH {
    tag "$meta.id"
    label 'process_single'

    container "labsyspharm/unetcoreograph:2.2.9"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*[0-9]*.tif"), emit: cores
    tuple val(meta), path("./masks/"), emit: masks
    tuple val(meta), path("TMA_MAP.tif"), emit: TMA_map
    tuple val(meta), path("centroidsY-X.txt"), emit: centroids

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /app/UNetCoreograph.py \\
        --imagePath ${image}\\
        --outputPath .\\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreograph: \$(echo 2.2.9)
    END_VERSIONS
    """
}
