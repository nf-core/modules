process DEEPTOOLS_COMPUTEMATRIXOPERATIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.mat.gz"), emit: matrix
    tuple val("${task.process}"), val('deeptools'), eval('deeptools --version | sed "s/deeptools //g"'), topic: versions, emit: versions_deeptools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''   // controls the main command - sort, relabel etc
    def args2  = task.ext.args2 ?: '' // accessory commands, needed to fine tune tool function 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    computeMatrixOperations \\
        ${args} \\
        -m ${matrix} \\
        -o ${prefix}.mat.gz \\
        ${args2}
    """

    stub:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    echo ${args2}
    echo "" | gzip > ${prefix}.mat.gz
    """
}
