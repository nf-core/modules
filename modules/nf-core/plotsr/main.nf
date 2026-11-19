// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLOTSR {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db8f402265ea63f1a4212d71f58474b4a72cba3b17f94469e9fda39bd0439ca1/data'
        : 'community.wave.seqera.io/library/plotsr_matplotlib_python:5179eb65b5a167ec'}"

    input:
    tuple val(meta), path(syri)
    tuple val(meta2), path(fastas)
    tuple val(meta3), path(genomes)
    tuple val(meta4), path(bedpe)
    tuple val(meta5), path(markers)
    tuple val(meta6), path(tracks)
    tuple val(meta7), path(chrord)
    tuple val(meta8), path(chrname)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val('plotsr'), eval("plotsr --version"), emit: versions_plotsr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def syri_arg = syri instanceof List ? syri.collect { syri_file -> "--sr ${syri_file}" }.join(' ') : "--sr ${syri}"
    def bedpe_arg = bedpe ? "--bedpe ${bedpe}" : ''
    def markers_arg = markers ? "--markers ${markers}" : ''
    def tracks_arg = tracks ? "--tracks ${tracks}" : ''
    def chrord_arg = chrord ? "--chrord ${chrord}" : ''
    def chrname_arg = chrname ? "--chrname ${chrname}" : ''
    """
    plotsr \\
        ${syri_arg} \\
        --genomes ${genomes} \\
        ${bedpe_arg} \\
        ${markers_arg} \\
        ${tracks_arg} \\
        ${chrord_arg} \\
        ${chrname_arg} \\
        ${args} \\
        -o ${prefix}.png
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    """
}
