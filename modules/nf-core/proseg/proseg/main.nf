process PROSEG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0dcce070a1e7b921edd0254596eb6945f97c54e2be0fe3130e2d2678b3cfd42/data':
        'community.wave.seqera.io/library/rust-proseg:2.0.4--6c02254be033edab' }"

    input:
    tuple val(meta), path(transcripts)
    val(mode)

    output: 
    tuple val(meta), path("transcript-metadata.csv.gz")     , emit: transcript_metadata
    tuple val(meta), path("cell-metadata.csv.gz")           , emit: cell_metadata
    tuple val(meta), path("cell-polygons.geojson.gz")       , emit: cell_polygons
    tuple val(meta), path("cell-polygons-layers.geojson.gz"), emit: cell_polygons_layers
    tuple val(meta), path("expected-counts.csv.gz")         , emit: expected_counts
    tuple val(meta), path("union-cell-polygons.geojson.gz") , emit: union_cell_polygons
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset = mode ? "--${mode}" : ''
    def args = task.ext.args ?: ''

    """
    proseg \\
        ${preset} \\
        ${args} \\
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
    touch union-cell-polygons.geojson.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed "s/proseg //")
    END_VERSIONS
    """
}
