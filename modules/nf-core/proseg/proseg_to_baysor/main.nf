process PROSEG_TO_BAYSOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38ebf0dd1e071eb5a99fd220459c09625c1465c5491a3e2dab392bfbe8acb45f/data':
        'community.wave.seqera.io/library/rust-proseg:2.0.5--dde937bdc1cf4715' }"

    input:
    tuple val(meta), path(transcript_metadata)
    tuple val(meta2), path(cell_polygons)

    output:
    tuple val(meta), path("*baysor-cell-polygons.geojson")  , emit: baysor_cell_polygons
    tuple val(meta), path("*baysor-transcript-metadata.csv"), emit: baysor_transcript_metadata
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    proseg-to-baysor \\
        ${transcript_metadata} \\
        ${cell_polygons} \\
        --output-transcript-metadata ${prefix}-baysor-transcript-metadata.csv \\
        --output-cell-polygons ${prefix}-baysor-cell-polygons.geojson \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-baysor-transcript-metadata.csv
    touch ${prefix}-baysor-cell-polygons.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """
}
