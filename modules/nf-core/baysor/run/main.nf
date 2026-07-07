process BAYSOR_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f716efdfa817ee57bf59b78248f672b503807bfb8608ad3d3976f2e706ca9fb4/data' :
        'community.wave.seqera.io/library/baysor:0.7.1--6fd896e03359bae6'}"

    input:
    tuple val(meta), path(transcripts), path(prior_segmentation), path(config), val(scale)
    val(prior_confidence)
    val(prior_column)
    val(polygon_format)

    output:
    tuple val(meta), path("${prefix}_segmentation.csv"), path("${prefix}_segmentation_polygons_2d.json"), emit: segmentation
    tuple val("${task.process}"), val('baysor'), eval("baysor --version"), topic: versions, emit: versions_baysor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def prior_seg = prior_column ? ":${prior_column}" : (prior_segmentation ?: '')
    def confidence = prior_confidence ? "--prior-segmentation-confidence=${prior_confidence}" : ''
    def scaling_factor = scale ? "--scale=${scale}" : ''
    def config_arg = config ? "--config=${config}" : ''

    // check for valid output polygon format
    def valid_formats = ['GeometryCollectionLegacy', 'GeometryCollection', 'FeatureCollection']
    if (!polygon_format in valid_formats) {
        error "Invalid output polygon format. Valid options: ${valid_formats.join(', ')}"
    }
    def polygon_fmt = polygon_format ? "--polygon-format=${polygon_format}" : '--polygon-format=FeatureCollection'

    """
    export JULIA_NUM_THREADS=${task.cpus}

    baysor \\
        run \\
        ${transcripts} \\
        ${prior_seg} \\
        ${scaling_factor} \\
        ${confidence} \\
        --output ${prefix}_segmentation.csv \\
        ${config_arg} \\
        ${polygon_fmt} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_segmentation.csv"
    touch "${prefix}_segmentation_polygons_2d.json"
    """
}
