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

    output: 
    tuple val(meta), path("*transcript-metadata.csv.gz", arity: '1')     , emit: transcript_metadata
    tuple val(meta), path("*cell-polygons-union.geojson.gz", arity: '1') , emit: union_cell_polygons 
    tuple val(meta), path("*cell-polygons.geojson.gz", arity: '1')       , emit: cell_polygons
    tuple val(meta), path("*cell-metadata.csv.gz", arity: '1')           , emit: cell_metadata
    tuple val(meta), path("*cell-polygons-layers.geojson.gz", arity: '1'), emit: cell_polygons_layers
    tuple val(meta), path("*expected-counts.csv.gz", arity: '1')         , emit: expected_counts
    tuple val(meta), path("*maxpost_counts*", arity: '1')                , emit: maxpost_counts, optional: true
    tuple val(meta), path("*output_rates*", arity: '1')                  , emit: output_rates, optional: true
    tuple val(meta), path("*cell-hulls*", arity: '1')                    , emit: cell_hulls, optional: true
    tuple val(meta), path("*gene-metadata*", arity: '1')                 , emit: gene_metadata, optional: true
    tuple val(meta), path("*metagene-rates*", arity: '1')                , emit: metagene_rates, optional: true
    tuple val(meta), path("*metagene-loadings*", arity: '1')             , emit: metagene_loadings, optional: true
    tuple val(meta), path("*cell-voxels*", arity: '1')                   , emit: cell_voxels, optional: true
    path "versions.yml"                                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def preset = mode ? "--${mode}" : ''
    def args = task.ext.args ?: ''
    def prefix =  task.ext.prefix ? "--output-path ${task.ext.prefix}" : ""

    """
    proseg \\
        ${preset} \\
        ${args} \\
        --output-union-cell-polygons cell-polygons-union.geojson.gz \\
        ${prefix} \\
        --nthreads ${task.cpus} \\
        ${transcripts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """

    stub:
    """
    touch cell-metadata.csv.gz
    touch cell-polygons.geojson.gz
    touch cell-polygons-layers.geojson.gz
    touch expected-counts.csv.gz
    touch transcript-metadata.csv.gz
    touch cell-polygons-union.geojson.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """
}
