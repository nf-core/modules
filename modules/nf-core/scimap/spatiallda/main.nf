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
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /scimap/scripts/spatialLDA.py \
        --input $phenotyped \
        --output "${prefix}.csv" \
        --neighborhood-composition-plot "${prefix}.png" \
        --motif-locations-plot "${prefix}.html" \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scimap/spatialLDA: \$(python /scimap/scripts/spatialLDA.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch "${prefix}.csv"
    touch "${prefix}.png"
    touch "${prefix}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scimap/spatialLDA: \$(python /scimap/scripts/spatialLDA.py --version)
    END_VERSIONS
    """
}
