// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FQ_LINT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1edbe701f591caa7f7a930781d0364392706d291016f1028b5b9fe501ec6ef8c/data'
        : 'community.wave.seqera.io/library/fq:0.12.0--b3082ff93d0ccff0'}"

    input:
    tuple val(meta), path(fastq, arity: '1..2')

    output:
    tuple val(meta), path("*.fq_lint.txt"), emit: lint
    tuple val("${task.process}"), val('fq'), eval("fq lint --version | sed 's/fq-lint //; s/ .*//'"), emit: versions_fq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq lint \\
        $args \\
        $fastq > ${prefix}.fq_lint.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fq_lint.txt
    """
}
