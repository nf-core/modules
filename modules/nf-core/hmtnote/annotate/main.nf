process HMTNOTE_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmtnote:0.7.2--pyhdfd78af_1':
        'quay.io/biocontainers/hmtnote:0.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_annotated.vcf"), emit: vcf
    tuple val("${task.process}"), val('hmtnote'), eval("hmtnote --version 2>&1 | sed 's/.*version //'"), emit: versions_hmtnote, topic: versions


    script:
    def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This tool is no longer maintained by its author and, as its database hosting service has
been discontinued, it can no longer work with conda.
"""
    assert false: deprecation_message
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    hmtnote \\
        annotate \\
        ${vcf} \\
        ${prefix}_annotated.vcf \\
        ${args}
    """

    stub:
        def deprecation_message = """
WARNING: This module has been deprecated.

Reason:
This tool is no longer maintained by its author and, as its database hosting service has
been discontinued, it can no longer work with conda.
"""
    assert false: deprecation_message
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.vcf
    """
}
