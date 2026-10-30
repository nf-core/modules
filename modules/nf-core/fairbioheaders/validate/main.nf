// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FAIRBIOHEADERS_VALIDATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d07f6d7a9ec156ba40d7b5c2b4f2f1685aea26601fdcec7252f3146fc61fceb0/data'
        : 'community.wave.seqera.io/library/fair-bioheaders:0.8.5--e659fa573f6f4d5b' }"

    input:
    tuple val(meta), path(metadata)

    output:
    tuple val(meta), path("${prefix}.validate.txt"), emit: report
    tuple val("${task.process}"), val('fair-bioheaders'), eval("bioheaders --version"), topic: versions, emit: versions_fairbioheaders

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    // A group, not a pipe, so that invalid metadata fails the task.
    """
    {
        printf 'file: %s\\n' ${quote.call(metadata)}
        bioheaders \\
            validate \\
            ${args} \\
            ${quote.call(metadata)}
    } > ${quote.call("${prefix}.validate.txt")}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def quote = { value -> "'" + value.toString().replace("'", "'\\''") + "'" }
    """
    touch ${quote.call("${prefix}.validate.txt")}
    """
}
