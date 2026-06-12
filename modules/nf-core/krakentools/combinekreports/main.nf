process KRAKENTOOLS_COMBINEKREPORTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2.1--pyh7e72e81_0':
        'quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(kreports)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('krakentools'), val('1.2.1'), topic: versions, emit: versions_krakentools


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    combine_kreports.py \\
        -r ${kreports} \\
        -o ${prefix}.txt \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
