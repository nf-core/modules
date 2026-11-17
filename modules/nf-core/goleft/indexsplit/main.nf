// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GOLEFT_INDEXSPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95b388cf83447fc133e87737248025aeca3e8f05bf545c4f5010ad553338ee0a/data'
        : 'community.wave.seqera.io/library/goleft:0.2.6--3aecf4b0548543d2'}"

    input:
    tuple val(meta) , path(bai)
    tuple val(meta2), path(fai)
    val split

    output:
    tuple val(meta), path("*.bed") , emit: bed
    tuple val("${task.process}"), val('goleft'), eval("goleft --version 2>&1 | sed '1!d;s/^.*goleft Version: //'"), topic: versions, emit: versions_goleft

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    goleft \\
        indexsplit \\
        ${args} \\
        -n ${split} \\
        --fai ${fai} \\
        ${bai} \\
        > ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
