process EPIC2_EPIC2 {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1d91a2e8076ff1f0f745df24a4603843e9af881673f6f4b06af2050242826918/data'
        : 'community.wave.seqera.io/library/epic2_python_setuptools:d5718ed0ebe5e13d'}"

    input:
    tuple val(meta), path(ipbam), path(controlbam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('epic2'), eval("epic2 --version"), topic: versions, emit: versions_epic2
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = controlbam ? "--control ${controlbam}" : ''
    """
    epic2 \\
        ${args} \\
        ${control} \\
        --treatment ${ipbam} \\
        -o ${prefix}.bed \\
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.bed
    """
}
