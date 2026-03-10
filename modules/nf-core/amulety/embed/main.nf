process AMULETY_EMBED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "docker.io/immcantation/amulety:2.1.2" // Seqera containers cannot be used since GPU is needed at runtime for pytorch with CUDA support to be installed

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
