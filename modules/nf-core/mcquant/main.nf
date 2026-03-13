process MCQUANT {
    tag "${meta.id}"
    label 'process_single'

    container "nf-core/quantification:1.5.4"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(mask)
    tuple val(meta3), path(markerfile)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('mcquant'), val("1.5.4"), emit: versions_mcquant, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python /app/CommandSingleCellExtraction.py \
        --masks ${mask} \
        --image ${image} \
        --channel_names ${markerfile} \
        --output . \
        ${args}
    """

    stub:
    """
    touch cycif_tonsil_registered_cell.csv
    """
}
