def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/dastool/fastatocontig2bin

Reason:
This tool has been renamed in newer versions of DAS_Tool, so any changes
to this tool will not be tracked by this module.
"""

process DASTOOL_SCAFFOLDS2BIN {
    tag "$meta.id"
    label 'process_single'

    // Do not bump! This is the 'old name' of contigs2bin which is only available up until 1.1.3!
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/das_tool:1.1.3--r41hdfd78af_0' :
        'biocontainers/das_tool:1.1.3--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val(extension)

    output:
    tuple val(meta), path("*.tsv"), emit: scaffolds2bin
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assert false: deprecation_message
    """
}
