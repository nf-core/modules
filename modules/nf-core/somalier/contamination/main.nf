// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SOMALIER_CONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f5/f5ad8a7223d537f18f2899b16bd2b14f93866656c4ddce6d16181319d9a82f63/data'
        : 'community.wave.seqera.io/library/somalier:0.3.6--127e0e29751f4940'}"

    input:
    tuple val(meta), path(extract, arity: '1..*')
    tuple val(meta2), path(sites)

    output:
    tuple val(meta), path("*.samples.tsv"), emit: samples_tsv
    tuple val(meta), path("*.pairs.tsv")  , emit: pairs_tsv
    tuple val("${task.process}"), val('somalier'), eval('somalier 2>&1 | sed -n \'s/.*version: \\([0-9.]*\\).*/\\1/p\''), emit: versions_somalier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    somalier contamination \\
        --sites ${sites} \\
        --output ${prefix} \\
        ${args} \\
        ${extract.join(' ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.samples.tsv
    touch ${prefix}.pairs.tsv
    """
}
