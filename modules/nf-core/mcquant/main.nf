process MCQUANT {
    tag "$meta.id"
    label 'process_single'

    container "labsyspharm/quantification:1.5.4"

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

    """
    python /app/CommandSingleCellExtraction.py \
        --masks $mask \
        --image $image \
        --output . \
        --channel_names $markerfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: \$(echo 1.5.4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch cycif_tonsil_registered_cell.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcquant: \$(echo 1.5.4)
    END_VERSIONS
    """
}
