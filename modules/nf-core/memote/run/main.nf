// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MEMOTE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6ccfddd278dd7116f7d952f138132b846da902d06fcf238a0a9e5595b632bf7/data'
        : 'community.wave.seqera.io/library/memote:0.17.0--206338ee32d58ea0'}"

    input:
    tuple val(meta), path(model)

    output:
    tuple val(meta), path("*.json.gz"), emit: json   , topic: json
    tuple val("${task.process}"), val('memote'), eval("memote --version | sed 's/memote, version //'"), topic: versions, emit: versions_memote

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\${PWD}
    export COBRA_SOLVER=glpk_exact

    memote run \\
        --filename ${prefix}.json.gz \\
        --ignore-git \\
        ${args} \\
        ${model}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\${PWD}
    export COBRA_SOLVER=glpk_exact

    echo "" | gzip > ${prefix}.json.gz
    """
}
