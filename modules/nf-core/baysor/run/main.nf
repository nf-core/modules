process BAYSOR_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e642d48a06693a7bad75c491ddf0d027ffa1450a718a9f3af573645662ebee48/data' :
        'community.wave.seqera.io/library/baysor_python:55a5c12eda9597fb' }"

    input:
    tuple val(meta), path(transcripts), path(prior_segmentation), path(config), val(scale)
    val(prior_confidence)
    val(prior_column)

    output:
    tuple val(meta), path("${prefix}/segmentation.csv"), path("${prefix}/segmentation_polygons_2d.json"), emit: segmentation
    tuple val("${task.process}"), val('baysor'), eval("baysor --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo unknown"), topic: versions, emit: versions_baysor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    // Column-based prior (e.g. :cell_id) takes precedence over file-based prior
    def prior_col = prior_column ? ":${prior_column}" : ''
    def prior_seg = prior_col ?: (prior_segmentation ? prior_segmentation : '')
    def confidence = prior_confidence != null ? "--prior-segmentation-confidence=${prior_confidence}" : ''
    def scaling_factor = scale ? "--scale=${scale}" : ''
    def config_arg = config ? "--config=${config}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export JULIA_NUM_THREADS=${task.cpus}

    mkdir -p ${prefix}

    baysor \\
        run \\
        ${transcripts} \\
        ${prior_seg} \\
        ${scaling_factor} \\
        ${confidence} \\
        --output="${prefix}/segmentation.csv" \\
        ${config_arg} \\
        --plot \\
        --polygon-format=GeometryCollectionLegacy \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/segmentation.csv"
    touch "${prefix}/segmentation_polygons_2d.json"
    """
}