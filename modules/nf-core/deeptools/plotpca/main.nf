// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEEPTOOLS_PLOTPCA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba9641d61acf134c6b989b3af5f25bcd5aadb3faf692a41c00ba4cd25f99032e/data'
:         'community.wave.seqera.io/library/deeptools:3.5.6--511a2b05870b9ec3' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('deeptools'), eval('plotPCA --version | sed "s/plotPCA //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotPCA \\
        $args \\
        --corData $matrix \\
        --plotFile ${prefix}.plotPCA.pdf \\
        --outFileNameData ${prefix}.plotPCA.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotPCA.pdf
    touch ${prefix}.plotPCA.tab
    """
}
