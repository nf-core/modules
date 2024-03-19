process SCIMAP_SPATIALLDA {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/miguelib/scimap:latest"

    input:
    tuple val(meta), path(phenotyped)

    output:
    tuple val(meta), path("*.csv"), emit: spatial_lda_output
    tuple val(meta), path("*.png"), emit: composition_plot
    tuple val(meta), path("*.html"), emit: motif_location_plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /scimap/scripts/spatialLDA.py \
        --input $phenotyped \
        --output $spatial_lda_output \
        --neighborhood-composition-plot $composition_plot \
        --motif-locations-plot $motif_location_plot \
        $args"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scimap/spatialLDA: \$(python /scimap/scripts/spatialLDA.py --version)
    END_VERSIONS
    """
}
