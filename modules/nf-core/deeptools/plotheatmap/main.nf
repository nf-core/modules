process DEEPTOOLS_PLOTHEATMAP {
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
    tuple val(meta), path("*.tab"), emit: table
    tuple val("${task.process}"), val('deeptools'), eval('plotHeatmap --version | sed "s/plotHeatmap //g"') , emit: versions_deeptools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotHeatmap \\
        $args \\
        --matrixFile $matrix \\
        --outFileName ${prefix}.plotHeatmap.pdf \\
        --outFileNameMatrix ${prefix}.plotHeatmap.mat.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotHeatmap.pdf
    touch ${prefix}.plotHeatmap.mat.tab
    """
}
