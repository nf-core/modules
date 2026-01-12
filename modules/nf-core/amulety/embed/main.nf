process AMULETY_EMBED {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/amulety_wget:57c081c80e888ece':
        'community.wave.seqera.io/library/amulety_wget:ff98ce74f1246e24' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)
    val(model)

    output:
    tuple val(meta), path("*.tsv"), emit: embedding
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TRANSFORMERS_CACHE="./cache" amulety \\
        embed \\
        $args \\
        --input-airr $tsv \\
        --chain $chain \\
        --model $model \\
        --output-file-path ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
    END_VERSIONS
    """
}
