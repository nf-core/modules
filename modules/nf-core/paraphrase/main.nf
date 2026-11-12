// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PARAPHRASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/067ee19a995c4594495da4a1098114e1b9e5390af80479b15d6ae68e67e5c6e1/data'
:         'community.wave.seqera.io/library/pip_python_paraphrase:48432dd9875d118c' }"

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
