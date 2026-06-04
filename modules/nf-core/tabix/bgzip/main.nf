process TABIX_BGZIP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("output.gz"), emit: output
    tuple val(meta), path("*.gzi")    , emit: gzi, optional: true
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")   , topic: versions   , emit: versions_tabix


    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use HTSLIB/BGZIPTABIX instead.

    Reason:
    This module is duplicative of TABIX/BGZIPTABIX and SAMTOOLS/BGZIP. The new HTSLIB/BGZIPTABIX module provides equivalent functionality with a more predictable behavior and better interface.
    """.stripIndent()
    assert false: deprecation_message

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use HTSLIB/BGZIPTABIX instead.

    Reason:
    This module is duplicative of TABIX/BGZIPTABIX and SAMTOOLS/BGZIP. The new HTSLIB/BGZIPTABIX module provides equivalent functionality with a more predictable behavior and better interface.
    """.stripIndent()
    assert false: deprecation_message
}
