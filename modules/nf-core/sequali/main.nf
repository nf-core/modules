process SEQUALI {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9ed3be5b0d3beb64807ec93b25a80b55abdcdffe684114d12ddef78461dd64e9/data':
        'community.wave.seqera.io/library/sequali:0.12.0--07485bec824d914a' }"

    input:

    tuple val(meta), path(reads)

    output:

    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_1_bam = reads.size() == 1 ? reads : reads[0]
    def read_2 = reads.size() == 2 ? reads[1]: ""

    """
    sequali \\
        $args \\
        -t $task.cpus \\
        --html ${prefix}.html \\
        --json ${prefix}.json \\
        $read_1_bam \\
        $read_2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequali: \$(sequali --version |& sed '1!d ; s/sequali //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.html
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequali: \$(sequali --version |& sed '1!d ; s/sequali //')
    END_VERSIONS
    """
}

