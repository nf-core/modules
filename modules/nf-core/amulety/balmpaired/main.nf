process AMULETY_BALMPAIRED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/amulety_wget:d69a2bc09a42a8b2':
        'community.wave.seqera.io/library/amulety_wget:2ecd2554d8d6f58e' }"

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
