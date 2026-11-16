// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_QUICKCHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/5749b5efca036c73d7b427c0ba347c94773330fe228b8eb395cdb619e7a5635c/data'
        : 'community.wave.seqera.io/library/samtools_htslib:1.25--65a8db39c6e3f5bd'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), env("EXIT_CODE"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools quickcheck \\
        ${args} \\
        ${bam} || EXIT_CODE=\$?

    EXIT_CODE=\${EXIT_CODE:-0}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    EXIT_CODE=0
    echo ${args}
    """
}
