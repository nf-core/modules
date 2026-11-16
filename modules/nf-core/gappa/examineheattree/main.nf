// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GAPPA_EXAMINEHEATTREE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4a30e4444b69c99470cae4a631a9e293ed76d475a6a44812a1f38ad2b9ece7af/data'
        : 'community.wave.seqera.io/library/gappa:0.9.0--2f0efa59bcb34481'}"

    input:
    tuple val(meta), path(jplace)

    output:
    tuple val(meta), path("*.newick")     , emit: newick  , optional: true
    tuple val(meta), path("*.nexus")      , emit: nexus   , optional: true
    tuple val(meta), path("*.phyloxml")   , emit: phyloxml, optional: true
    tuple val(meta), path("*.svg")        , emit: svg     , optional: true
    tuple val(meta), path("*.colours.txt"), emit: colours
    tuple val(meta), path("*.log")        , emit: log
    tuple val("${task.process}"), val('gappa'), eval("gappa --version 2>&1 | sed 's/v//'"), emit: versions_gappa, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gappa \\
        examine \\
        heat-tree \\
        --threads ${task.cpus} \\
        --file-prefix ${prefix}. \\
        --jplace-path ${jplace} \\
        ${args} \\
        --log-file ${prefix}.log

    grep '^ *At' ${prefix}.log > ${prefix}.colours.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.colours.txt
    touch ${prefix}.log
    """
}
