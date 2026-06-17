process SAMTOOLS_BGZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(infile)
    val out_ext

    output:
    tuple val(meta), path("output.gz"), emit: output
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions


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
