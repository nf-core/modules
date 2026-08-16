process BUSCO_PLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/963bad66c10646cf0adb1967cc462ad04d02789ddbfae4fbb94182291dbddf8c/data'
        : 'community.wave.seqera.io/library/busco:6.1.0--6d1f7006d91892b3'}"

    input:
    path short_summary_json, stageAs: 'busco/*'

    output:
    path '*.png', emit: png
    tuple val("${task.process}"), val('busco'), eval("busco --version | sed 's/BUSCO //'"), emit: versions_busco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'busco_figure'
    """
    busco \\
        --plot busco \\
        ${args}

    mv ./busco/busco_figure.png ${prefix}.png
    """

    stub:
    def prefix = task.ext.prefix ?: 'busco_figure'
    """
    touch ${prefix}.png
    """
}
