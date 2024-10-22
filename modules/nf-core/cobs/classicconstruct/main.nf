process COBS_CLASSICCONSTRUCT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cobs:0.3.0--hdcf5f25_1' :
        'biocontainers/cobs:0.3.0--hdcf5f25_1'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.index.cobs_classic")   , emit: index
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def task_memory_in_bytes = task.memory.toBytes()
    """
    cobs \\
        classic-construct \\
        $args \\
        --memory $task_memory_in_bytes \\
        --threads $task.cpus \\
        $input \\
        ${prefix}.index.cobs_classic


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobs: \$(cobs version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.index.cobs_classic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobs: \$(cobs version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
