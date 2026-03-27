process GLIMPSE_CHUNK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--h2ce4488_2'
        : 'biocontainers/glimpse-bio:1.1.1--hce55b13_1'}"

    input:
    tuple val(meta), path(input), path(input_index), val(region)

    output:
    tuple val(meta), path("*.txt"), emit: chunk_chr
    tuple val("${task.process}"), val('glimpse'), eval("GLIMPSE_chunk --help | sed -n '/Version/s/.*: //p'"), topic: versions, emit: versions_glimpse

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""

    """
    GLIMPSE_chunk \\
        ${args} \\
        --input ${input} \\
        --region ${region} \\
        --thread ${task.cpus} \\
        --output ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    """
}
