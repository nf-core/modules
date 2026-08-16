process SCIMAP_SPATIALLDA {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/miguelib/scimap:0.0.3"

    input:
    tuple val(meta), path(phenotyped)

    output:
    tuple val(meta), path("*.csv") , emit: spatial_lda_output
    tuple val(meta), path("*.png") , emit: composition_plot
    tuple val(meta), path("*.html"), emit: motif_location_plot
    tuple val("${task.process}"), val('scimap'), eval('python /scimap/scripts/spatialLDA.py --version'), emit: versions_scimap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${phenotyped}" == "${prefix}.csv") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate."

    """
    python /scimap/scripts/spatialLDA.py \\
        --input $phenotyped \\
        --output "${prefix}.csv" \\
        --neighborhood-composition-plot "${prefix}.png" \\
        --motif-locations-plot "${prefix}.html" \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch "${prefix}.csv"
    touch "${prefix}.png"
    touch "${prefix}.html"
    """
}
