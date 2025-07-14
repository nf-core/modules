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
    tuple val(meta), path("outputs/*transcript-metadata.csv.gz")     , emit: transcript_metadata
    tuple val(meta), path("outputs/*cell-metadata.csv.gz")           , emit: cell_metadata
    tuple val(meta), path("outputs/*cell-polygons.geojson.gz")       , emit: cell_polygons
    tuple val(meta), path("outputs/*cell-polygons-layers.geojson.gz"), emit: cell_polygons_layers
    tuple val(meta), path("outputs/*expected-counts.csv.gz")         , emit: expected_counts
    tuple val(meta), path("outputs/*")                               , emit: other_outputs, optional: true
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset = mode ? "--${mode}" : ''
    def args = task.ext.args ?: ''
    def prefix =  task.ext.prefix ? "outputs/${task.ext.prefix}-" : "outputs/"

    """
    mkdir outputs
    proseg \\
        ${preset} \\
        ${args} \\
        --output-transcript-metadata ${prefix}transcript-metadata.csv.gz \\
        --output-cell-metadata ${prefix}cell-metadata.csv.gz \\
        --output-cell-polygons ${prefix}cell-polygons.geojson.gz \\
        --output-cell-polygon-layers ${prefix}cell-polygons-layers.geojson.gz \\
        --output-expected-counts ${prefix}expected-counts.csv.gz \\
        --nthreads ${task.cpus} \\
        ${transcripts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """

    stub:
    """
    mkdir outputs
    touch outputs/cell-metadata.csv.gz
    touch outputs/cell-polygons.geojson.gz
    touch outputs/cell-polygons-layers.geojson.gz
    touch outputs/expected-counts.csv.gz
    touch outputs/transcript-metadata.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """
}
