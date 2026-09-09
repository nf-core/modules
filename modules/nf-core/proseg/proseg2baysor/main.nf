process PROSEG2BAYSOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c9a638dd3ab6d2ce5f366a8aaec27dc34dcf2914f01112a22b855c71bf58a7b/data' :
        'community.wave.seqera.io/library/rust-proseg:3.1.1--5dbb81bc6361cb10'}"

    input:
    tuple val(meta), path(sd_zarr)

    output:
    tuple val(meta), path("${prefix}_cell-polygons.geojson")  , emit: cell_polygons
    tuple val(meta), path("${prefix}_transcript-metadata.csv"), emit: transcript_metadata
    tuple val("${task.process}"), val('proseg'), eval("proseg --version | sed 's/proseg //'"), topic: versions, emit: versions_proseg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    proseg-to-baysor \\
        --output-transcript-metadata ${prefix}_transcript-metadata.csv \\
        --output-cell-polygons ${prefix}_cell-polygons.geojson \\
        ${sd_zarr} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_transcript-metadata.csv
    touch ${prefix}_cell-polygons.geojson
    """
}
