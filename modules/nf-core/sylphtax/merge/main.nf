// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SYLPHTAX_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f299731866f7b830801a515c1f624e6fabf5fe2e22af31d347409636c21db3b9/data'
        : 'community.wave.seqera.io/library/sylph-tax:1.9.2--861d9f8322018f98'}"

    input:
    tuple val(meta), path(sylphtax_reports)
    val data_type

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('sylph-tax'), eval("sylph-tax --version 2>&1 | tail -1"), emit: versions_sylphtax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    sylph-tax \\
        merge \\
        ${sylphtax_reports} \\
        --column ${data_type} \\
        --output ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    touch ${prefix}.tsv
    """
}
