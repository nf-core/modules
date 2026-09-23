process AMULETY_ESM2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/ae2be6d5fd5a4e1024b78eb1acdfdcb6aab4326f002e98d7c2c97ef00aa979e2/data'
:         'community.wave.seqera.io/library/amulety:1.1--5abbe5fc5e136fd3' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)

    output:
    tuple val(meta), path("*.tsv"), emit: embedding
    tuple val("${task.process}"), val('amulety'), eval("amulety --help 2>&1 | grep -o 'version [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_amulety, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/amulety/embed instead

Reason:
This module is no longer fit for purpose because the syntax for amulety has been updated in version 2.x.
The new 'embed' command now covers the embedding functionality for all embeddings.

"""
    assert false: deprecation_message

    stub:
    def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/amulety/embed instead

Reason:
This module is no longer fit for purpose because the syntax for amulety has been updated in version 2.x.
The new 'embed' command now covers the embedding functionality for all embeddings.

"""
    assert false: deprecation_message
}
