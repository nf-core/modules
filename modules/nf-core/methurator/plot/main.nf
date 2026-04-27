process METHURATOR_PLOT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/methurator:2.1.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/methurator:2.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(summary_report)

    output:
    tuple val(meta), path("plots/*.html"), emit: plots
    tuple val("${task.process}"), val('methurator'), eval("methurator --version | sed 's/.* //'"), emit: versions_methurator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    methurator plot \\
        --summary ${summary_report} \\
        --outdir .

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir plots/
    touch plots/${prefix}.html

    """
}
