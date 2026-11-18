// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KAIJU_MERGEOUTPUTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/40cd79a6bb3b5feabe877bef8498d63a3b91fed44b30556af4c7e6554c4310f2/data'
        : 'community.wave.seqera.io/library/kaiju:1.10.3--8b4a2069cc29ac7a'}"

    input:
    tuple val(meta), path(kaiju), path(kraken)
    path db

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    tuple val("${task.process}"), val('kaiju'), eval("kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'"), emit: versions_kaiju, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbnodes = db ? '-t <(find -L ${db} -name "*nodes.dmp")' : ''

    if ("${kaiju}" == "${prefix}.tsv") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    if ("${kraken}" == "${prefix}.tsv") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    kaiju-mergeOutputs \\
        -i <(sort -k2,2 ${kaiju}) \\
        -j <(sort -k2,2 ${kraken}) \\
        -o ${prefix}.tsv \\
        ${dbnodes} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
