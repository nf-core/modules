// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_BGZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4'}"

    input:
    tuple val(meta), path(infile)
    val out_ext

    output:
    tuple val(meta), path("output.gz"), emit: output
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use HTSLIB/BGZIPTABIX instead.

    Reason:
    This module is duplicative of TABIX/BGZIPTABIX and HTSLIB/BGZIPTABIX. The new HTSLIB/BGZIPTABIX module provides equivalent functionality with a more predictable behavior and better interface.
    """.stripIndent()
    assert false: deprecation_message

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use HTSLIB/BGZIPTABIX instead.

    Reason:
    This module is duplicative of TABIX/BGZIPTABIX and HTSLIB/BGZIPTABIX. The new HTSLIB/BGZIPTABIX module provides equivalent functionality with a more predictable behavior and better interface.
    """.stripIndent()
    assert false: deprecation_message
}
