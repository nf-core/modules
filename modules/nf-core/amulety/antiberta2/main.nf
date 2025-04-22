process AMULETY_ANTIBERTA2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/amulety_igblast:b2a7736f645c40e5':
        'community.wave.seqera.io/library/amulety_igblast:659eaa872785adeb' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)

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
        antiberta2 \\
        $args \\
        --cache-dir ./cache \\
        $tsv \\
        $chain \\
        ${prefix}.tsv

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
