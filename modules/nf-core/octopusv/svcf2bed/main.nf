process OCTOPUSV_SVCF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8378c503810508764badf4504647ea1791fa1f585dc5b77ad5710237609b15ff/data'
:         'community.wave.seqera.io/library/octopusv:1.0.0--5a8c89c65098e801' }"

    input:
    tuple val(meta), path(svcf)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    octopusv svcf2bed \
        --input-file ${svcf} \\
        --output-file ${prefix}.bed \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.bed
    """
}
