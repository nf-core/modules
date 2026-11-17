// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SPATIALSANCTUARY_SIMULATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2d3e7b96206df4ad45ef3374ae0c54d27077474fc78502cbfd0f8bca684f7504/data'
        : 'community.wave.seqera.io/library/spatial-sanctuary:2.2.2--da11ed94efb15111'}"

    input:
    tuple val(meta), val(modality), val(c0), val(theta), val(r_ox)

    output:
    tuple val(meta), path("*_sanctuary.json"), emit: json
    tuple val("${task.process}"), val('spatial-sanctuary'), eval("python3 -c \"import importlib.metadata; print(importlib.metadata.version('spatial-sanctuary'))\""), topic: versions, emit: versions_spatialsanctuary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    spatial-sanctuary simulate \\
        --modality "${modality}" \\
        --c0 "${c0}" \\
        --theta "${theta}" \\
        --r-ox "${r_ox}" \\
        --json \\
        --out "${prefix}_sanctuary.json" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_sanctuary.json"
    """
}
