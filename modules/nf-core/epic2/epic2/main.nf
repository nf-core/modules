process EPIC2_EPIC2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/22722cc6677e4fdec92038f59d73081b77e18e040bc12870eb56c49844a9c5f8/data':
        'community.wave.seqera.io/library/epic2_python:869ce1d0bdd993e2' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('epic2'), eval("epic2 --version"), topic: versions, emit: versions_epic2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = controlbam ? "--control $controlbam" : ''

    """
    epic2 \\
        $args \\
        $control \\
        --treatment $ipbam \\
        -o ${prefix}.bed \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.bed
    """
}
