// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CADDSV_GET {
    tag "CADDSV ${flag}"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1f3e808e9aa641e613532f33d93d0b3f6b242a7059d5ededc9cd9df0d4d7b2fd/data'
        : 'community.wave.seqera.io/library/caddsv:2.0.4--888cae10b93818a7'}"

    input:
    val flag

    output:
    path "caddsv_annotations", emit: annotations
    tuple val("${task.process}"), val('caddsv'), eval("caddsv --version"), emit: versions_caddsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (!['annotations', 'segmentnt'].contains(flag)) {
        error("Invalid caddsv get flag: '${flag}'. Expected 'annotations' or 'segmentnt'.")
    }
    """
    caddsv get ${flag} \\
        --annotations-dir "caddsv_annotations" \\
        ${args}
    """

    stub:
    """
    mkdir -p caddsv_annotations
    touch caddsv_annotations/stub.txt
    """
}
