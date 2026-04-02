process AMULETY_EMBED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/eef2fbc7c8d1ba71b3890a83b520c3eefa136ec4de5a8e6a97db828ae354d7ab/data':
        'community.wave.seqera.io/library/igblast_curl_python_wget_pruned:07dda71433b05ed5' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)
    val(model)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.tsv"), emit: embedding
    tuple val(meta), path("*metadata.tsv"), emit: embedding_metadata
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
    touch ${prefix}_metadata.tsv
    """
}
