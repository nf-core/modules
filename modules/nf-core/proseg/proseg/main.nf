process PROSEG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c9a638dd3ab6d2ce5f366a8aaec27dc34dcf2914f01112a22b855c71bf58a7b/data' :
        'community.wave.seqera.io/library/rust-proseg:3.1.1--5dbb81bc6361cb10'}"

    input:
    tuple val(meta), path(transcripts)
    val mode
    tuple val(transcript_metadata_fmt), val(cell_metadata_fmt), val(expected_counts_fmt)

    output:
    tuple val(meta), path("${prefix}.zarr"                                    , arity: '1'), emit: zarr
    tuple val(meta), path("${prefix}_transcript-metadata.{csv,csv.gz,parquet}", arity: '1'), emit: transcript_metadata
    tuple val(meta), path("${prefix}_cell-polygons.geojson.gz"                , arity: '1'), emit: cell_polygons
    tuple val(meta), path("${prefix}_cell-metadata.{csv,csv.gz,parquet}"      , arity: '1'), emit: cell_metadata
    tuple val(meta), path("${prefix}_cell-polygons-layers.geojson.gz"         , arity: '1'), emit: cell_polygons_layers
    tuple val(meta), path("${prefix}_expected-counts.{csv,csv.gz,parquet}"    , arity: '1'), emit: expected_counts
    tuple val(meta), path("${prefix}_cell-polygons-union.geojson.gz"          , arity: '1'), emit: union_cell_polygons
    tuple val(meta), path("${prefix}_maxpost-counts.{csv,csv.gz,parquet}"                 ), emit: maxpost_counts      , optional: true
    tuple val(meta), path("${prefix}_output-rates.{csv,csv.gz,parquet}"                   ), emit: output_rates        , optional: true
    tuple val(meta), path("${prefix}_cell-hulls.{csv,csv.gz,parquet}"                     ), emit: cell_hulls          , optional: true
    tuple val(meta), path("${prefix}_gene-metadata.{csv,csv.gz,parquet}"                  ), emit: gene_metadata       , optional: true
    tuple val(meta), path("${prefix}_metagene-rates.{csv,csv.gz,parquet}"                 ), emit: metagene_rates      , optional: true
    tuple val(meta), path("${prefix}_metagene-loadings.{csv,csv.gz,parquet}"              ), emit: metagene_loadings   , optional: true
    tuple val(meta), path("${prefix}_cell-voxels.{csv,csv.gz,parquet}"                    ), emit: cell_voxels         , optional: true
    tuple val("${task.process}"), val('proseg'), eval("proseg --version | sed 's/proseg //'"), topic: versions, emit: versions_proseg


    script:
    def preset = mode ? "--${mode}" : ''
    def args = task.ext.args ?: ''
    prefix =  task.ext.prefix ?: "${meta.id}"

    def trans_meta_fmt = transcript_metadata_fmt ?: 'csv.gz'
    def cell_meta_fmt = cell_metadata_fmt ?: 'csv.gz'
    def exp_counts_fmt = expected_counts_fmt ?: 'csv.gz'

    def output_formats = ['csv', 'csv-gz', 'parquet']
    if (!trans_meta_fmt in output_formats) {
        error "Unsupported transcript metadata output format: ${trans_meta_fmt}. Supported formats are: ${output_formats.join(', ')}"
    }
    if (!cell_meta_fmt in output_formats) {
        error "Unsupported cell metadata output format: ${cell_meta_fmt}. Supported formats are: ${output_formats.join(', ')}"
    }
    if (!exp_counts_fmt in output_formats) {
        error "Unsupported expected counts output format: ${exp_counts_fmt}. Supported formats are: ${output_formats.join(', ')}"
    }

    """
    proseg \\
        ${preset} \\
        ${args} \\
        --output-spatialdata ${prefix}.zarr \\
        --output-transcript-metadata ${prefix}_transcript-metadata.${trans_meta_fmt} \\
        --output-cell-polygons ${prefix}_cell-polygons.geojson.gz \\
        --output-cell-metadata ${prefix}_cell-metadata.${cell_meta_fmt} \\
        --output-expected-counts ${prefix}_expected-counts.${exp_counts_fmt} \\
        --output-cell-polygon-layers ${prefix}_cell-polygons-layers.geojson.gz \\
        --output-union-cell-polygons ${prefix}_cell-polygons-union.geojson.gz \\
        --nthreads ${task.cpus} \\
        ${transcripts}
    """

    stub:
    prefix =  task.ext.prefix ?: "${meta.id}"
    def trans_meta_fmt = transcript_metadata_fmt ?: 'csv.gz'
    def cell_meta_fmt = cell_metadata_fmt ?: 'csv.gz'
    def exp_counts_fmt = expected_counts_fmt ?: 'csv.gz'

    """
    touch ${prefix}.zarr
    echo "" | gzip > ${prefix}_cell-metadata.${cell_meta_fmt}
    echo "" | gzip > ${prefix}_cell-polygons.geojson.gz
    echo "" | gzip > ${prefix}_cell-polygons-layers.geojson.gz
    echo "" | gzip > ${prefix}_expected-counts.${exp_counts_fmt}
    echo "" | gzip > ${prefix}_transcript-metadata.${trans_meta_fmt}
    echo "" | gzip > ${prefix}_cell-polygons-union.geojson.gz
    """
}
