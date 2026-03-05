process PARAPHRASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc177b8e7e1d9bbdcbb64a6ad630c7ecc63a2229ce2b219408888bc2bb34cac3/data':
        'community.wave.seqera.io/library/pip_paraphrase:59a4576966ee5f0b' }"

    input:
    tuple val(meta), path(jsons), val(samples)
    tuple val(meta2), path(yaml)
    val(tsv_output)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true

    tuple val("${task.process}"), val('paraphrase'), eval("paraphrase --version | sed 's/.* //'"), topic: versions, emit: versions_paraphrase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rules = yaml ? "--rules $yaml" : ''
    def output_format = tsv_output ? 'tsv' : 'json'
    """
    paraphrase \
        $args \
        --input ${jsons.join(' --input ')} \
        --sample ${samples.join(' --sample ')} \
        --output-format=${output_format} \
        $rules \
        > ${prefix}.${output_format}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_format = tsv_output ? 'tsv' : 'json'
    """
    echo $args

    touch ${prefix}.${output_format}
    """
}
