process TIDK_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.7--h6872113_0':
        'quay.io/biocontainers/tidk:0.2.7--h6872113_0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    tuple val("${task.process}"), val('tidk'), eval("tidk --version | sed 's/tidk //'"), emit: versions_tidk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tidk \\
        plot \\
        --output $prefix \\
        $args \\
        --tsv "$tsv"

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svg

    """
}
