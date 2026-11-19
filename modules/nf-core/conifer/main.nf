// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CONIFER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/cc53c94256fb2b27b0b0a17c5aaa7629bb5f9fb75ba127cf6eaaba93943c9dd7/data'
        : 'community.wave.seqera.io/library/conifer:1.0.2--d9e7201487c31241'}"

    input:
    tuple val(meta), path(kraken_result)
    path kraken_taxon_db

    output:
    tuple val(meta), path("*.score"), emit: score
    tuple val("${task.process}"), val('conifer'), eval("conifer --version 2>&1 | sed 's/^.*Conifer //'"), topic: versions, emit: versions_conifer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    conifer \\
        $args \\
        --input $kraken_result \\
        --db $kraken_taxon_db > ${prefix}.score
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.score
    """
}
