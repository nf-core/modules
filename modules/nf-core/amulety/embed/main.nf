process AMULETY_EMBED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c39fc87288811f7806452ecbdb559b9e9bba71aebb82c60d60af939a73bdf614/data':
        'community.wave.seqera.io/library/igblast_curl_python_transformers_pruned:05685e2c81024d42' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)
    val(model)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.tsv"), emit: embedding
    tuple val(meta), path("*metadata.tsv"), emit: embedding_metadata
    tuple val("${task.process}"), val('amulety'), eval("amulety --help 2>&1 | grep -o 'version [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_amulety, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | grep -o 'Python [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pytorch'), eval("python -c 'import torch; print(torch.__version__)'"), emit: versions_pytorch, topic: versions
    tuple val("${task.process}"), val('transformers'), eval("python -c 'import transformers; print(transformers.__version__)'"), emit: versions_transformers, topic: versions

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
    touch ${prefix}_metadata.tsv
    """
}
