process SEQUALI {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequali:1.0.2--py312h0fa9677_0':
        'biocontainers/sequali:1.0.2--py312h0fa9677_0' }"

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
