process PROSEG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38ebf0dd1e071eb5a99fd220459c09625c1465c5491a3e2dab392bfbe8acb45f/data':
        'community.wave.seqera.io/library/rust-proseg:2.0.5--dde937bdc1cf4715' }"

    input:
    tuple val(meta), path(transcripts)
    val mode
    tuple val(transcript_metadata_fmt), val(cell_metadata_fmt), val(expected_counts_fmt)

    output:
    tuple val(meta), path("*transcript-metadata.{csv,csv.gz,parquet}", arity: '1'), emit: transcript_metadata
    tuple val(meta), path("*cell-polygons.geojson.gz"                , arity: '1'), emit: cell_polygons
    tuple val(meta), path("*cell-metadata.{csv,csv.gz,parquet}"      , arity: '1'), emit: cell_metadata
    tuple val(meta), path("*cell-polygons-layers.geojson.gz"         , arity: '1'), emit: cell_polygons_layers
    tuple val(meta), path("*expected-counts.{csv,csv.gz,parquet}"    , arity: '1'), emit: expected_counts
    tuple val(meta), path("*cell-polygons-union.geojson.gz"          , arity: '1'), emit: union_cell_polygons
    tuple val(meta), path("*maxpost-counts.{csv,csv.gz,parquet}"                 ), emit: maxpost_counts      , optional: true
    tuple val(meta), path("*output-rates.{csv,csv.gz,parquet}"                   ), emit: output_rates        , optional: true
    tuple val(meta), path("*cell-hulls.{csv,csv.gz,parquet}"                     ), emit: cell_hulls          , optional: true
    tuple val(meta), path("*gene-metadata.{csv,csv.gz,parquet}"                  ), emit: gene_metadata       , optional: true
    tuple val(meta), path("*metagene-rates.{csv,csv.gz,parquet}"                 ), emit: metagene_rates      , optional: true
    tuple val(meta), path("*metagene-loadings.{csv,csv.gz,parquet}"              ), emit: metagene_loadings   , optional: true
    tuple val(meta), path("*cell-voxels.{csv,csv.gz,parquet}"                    ), emit: cell_voxels         , optional: true
    path "versions.yml"                                                           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def preset = mode ? "--${mode}" : ''
    def args = task.ext.args ?: ''
    def prefix =  task.ext.prefix ?: "${meta.id}-"

    """
    proseg \\
        ${preset} \\
        ${args} \\
        --output-transcript-metadata ${prefix}transcript-metadata.${transcript_metadata_fmt} \\
        --output-cell-polygons ${prefix}cell-polygons.geojson.gz \\
        --output-cell-metadata ${prefix}cell-metadata.${cell_metadata_fmt} \\
        --output-expected-counts ${prefix}expected-counts.${expected_counts_fmt} \\
        --output-cell-polygon-layers ${prefix}cell-polygons-layers.geojson.gz \\
        --output-union-cell-polygons ${prefix}cell-polygons-union.geojson.gz \\
        --nthreads ${task.cpus} \\
        ${transcripts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """

    stub:
    def prefix =  task.ext.prefix ?: "${meta.id}-"
    """
    echo | gzip > ${prefix}cell-metadata.${cell_metadata_fmt}
    echo | gzip > ${prefix}cell-polygons.geojson.gz
    echo | gzip > ${prefix}cell-polygons-layers.geojson.gz
    echo | gzip > ${prefix}expected-counts.${expected_counts_fmt}
    echo | gzip > ${prefix}transcript-metadata.${transcript_metadata_fmt}
    echo | gzip > ${prefix}cell-polygons-union.geojson.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """
}
