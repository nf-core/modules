process MCQUANT {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "nf-core/quantification:1.5.4"

    input:
    tuple val(meta), path(image)
    tuple val(meta2), path(mask)
    tuple val(meta3), path(markerfile)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('mcquant'), eval("echo 1.5.4"), emit: versions_mcquant, topic: versions // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python /app/CommandSingleCellExtraction.py \
        --masks $mask \
        --image $image \
        --channel_names $markerfile \
        --output . \
        $args
    """

    stub:
    """
    touch cycif_tonsil_registered_cell.csv
    """
}
