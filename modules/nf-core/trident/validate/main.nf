process TRIDENT_VALIDATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/poseidon-trident:2.2.2.1--hf7d7819_0'
        : 'quay.io/biocontainers/poseidon-trident:2.2.2.1--hf7d7819_0'}"

    input:
    // Inputs need to be synchronised so that all files for a given package are provided together.
    tuple val(meta), path(input_package)

    output:
    tuple val("${task.process}"), val('trident'), eval('trident --version'), emit: versions_trident, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    trident validate \
        -d . \
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo trident validate -d . ${args}
    """
}
