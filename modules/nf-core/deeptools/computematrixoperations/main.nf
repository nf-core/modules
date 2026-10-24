// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEEPTOOLS_COMPUTEMATRIXOPERATIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba9641d61acf134c6b989b3af5f25bcd5aadb3faf692a41c00ba4cd25f99032e/data'
:         'community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3' }"

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
