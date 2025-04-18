process AMULETY_ANTIBERTA2 {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-92ebbfc09fc136b8e201cb187cd9567ba335d439:459e6ebe51fb2818cb6de807f2c5fa99599b1214-0':
        'biocontainers/mulled-v2-92ebbfc09fc136b8e201cb187cd9567ba335d439:459e6ebe51fb2818cb6de807f2c5fa99599b1214-0' }"

    input:
    tuple val(meta), path(tsv)
    val(chain)

    output:
    tuple val(meta), path("*_antiberta2.tsv"), emit: embedding
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TRANSFORMERS_CACHE="./cache"
    amulety \\
        antiberta2 \\
        $args \\
        --cache-dir ./cache \\
        $tsv \\
        $chain \\
        ${prefix}_antiberta2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_antiberta2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amulety: \$( amulety --help 2>&1 | grep -o "version [0-9\\.]\\+" | grep -o "[0-9\\.]\\+" )
    END_VERSIONS
    """
}
