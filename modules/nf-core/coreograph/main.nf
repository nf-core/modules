process COREOGRAPH {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/labsyspharm/unetcoreograph:2.2.9"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*[0-9]*.tif")     , emit: cores
    tuple val(meta), path("./masks/")        , emit: masks
    tuple val(meta), path("TMA_MAP.tif")     , emit: tma_map
    tuple val(meta), path("centroidsY-X.txt"), emit: centroids

    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    python /app/UNetCoreograph.py \\
        --imagePath ${image}\\
        --outputPath .\\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreograph: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '2.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch TMA_MAP.tif
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreograph: $VERSION
    END_VERSIONS
    """
}
