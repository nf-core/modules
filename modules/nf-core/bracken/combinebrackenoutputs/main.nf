process BRACKEN_COMBINEBRACKENOUTPUTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f30aa99d8d4f6ff1104f56dbacac95c1dc0905578fb250c80f145b6e80703bd1/data':
        'community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('combine_bracken_outputs'), eval('bracken -v | cut -f2 -d"v"'), emit: versions_combine_bracken_outputs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    combine_bracken_outputs.py \\
        $args \\
        --files ${input} \\
        -o ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
