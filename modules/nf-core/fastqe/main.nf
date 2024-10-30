process FASTQE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqe:0.3.3--pyhdfd78af_0':
        'biocontainers/fastqe:0.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqe \\
        $reads \\
        $args \\
        --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqe: \$(fastqe --version 2>&1 | sed 's/fastqe, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqe: \$(fastqe --version 2>&1 | sed 's/fastqe, version //')
    END_VERSIONS
    """
}
