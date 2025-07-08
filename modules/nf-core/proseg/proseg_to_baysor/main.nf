process PROSEG_TO_BAYSOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0dcce070a1e7b921edd0254596eb6945f97c54e2be0fe3130e2d2678b3cfd42/data':
        'community.wave.seqera.io/library/rust-proseg:2.0.4--6c02254be033edab' }"

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
    touch ${prefix}-transcript-metadata.csv
    touch ${prefix}-cell-polygons.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """
}
