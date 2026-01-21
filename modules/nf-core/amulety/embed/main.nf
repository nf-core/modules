process AMULETY_EMBED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/amulety_wget:6bee673d7a6a9753':
        'community.wave.seqera.io/library/amulety_wget:662e66d72e77dc3a' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)
    val(model)

    output:
    tuple val(meta), path("*.tsv"), emit: embedding
    tuple val("${task.process}"), val('amulety'), eval("amulety --help 2>&1 | grep -o 'version [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_amulety, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${tsv}" == "${prefix}.tsv") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    TRANSFORMERS_CACHE="./cache" amulety \\
        embed \\
        $args \\
        --input-airr $tsv \\
        --chain $chain \\
        --model $model \\
        --output-file-path ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${tsv}" == "${prefix}.tsv") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.tsv
    """
}
