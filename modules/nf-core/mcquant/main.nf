process MCQUANT {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "docker.io/labsyspharm/quantification:1.5.4"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(mask)
    tuple val(meta3), path(markerfile)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.5.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    python /app/CommandSingleCellExtraction.py \
        --masks $mask \
        --image $image \
        --channel_names $markerfile \
        --output . \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.5.4'
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: $VERSION
    END_VERSIONS
    """
}
