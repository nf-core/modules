// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTAGUARD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e76cd33dcff373007f86b4704c30c6aa0005b861debe63a9293fb15b577f9f8d/data'
        : 'community.wave.seqera.io/library/fastaguard:1.1.0--7f884e70d672a16c'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fastaguard.html"), emit: html
    tuple val(meta), path("*.fastaguard.json"), emit: json
    tuple val(meta), path("*.fastaguard.tsv"), emit: tsv
    tuple val(meta), path("*.fastaguard_mqc.json"), emit: mqc
    tuple val("${task.process}"), val('fastaguard'), eval('fastaguard --version | cut -d " " -f 2'), emit: versions_fastaguard, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    fastaguard ${fasta} \
      ${args} \
      --out ${prefix}.fastaguard.html \
      --json ${prefix}.fastaguard.json \
      --tsv ${prefix}.fastaguard.tsv \
      --multiqc ${prefix}.fastaguard_mqc.json
    """
}
