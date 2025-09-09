process NONPAREIL_NONPAREIL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.5.5--r43hdcf5f25_0':
        'biocontainers/nonpareil:3.5.5--r43hdcf5f25_0' }"

    input:
    tuple val(meta), path(reads)
    val format
    val mode

    output:
    tuple val(meta), path("*.npa"), emit: npa
    tuple val(meta), path("*.npc"), emit: npc
    tuple val(meta), path("*.npl"), emit: npl
    tuple val(meta), path("*.npo"), emit: npo
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    nonpareil \\
        -s $reads \\
        -f $format \\
        -T ${mode} \\
        -t $task.cpus \\
        -R ${mem_mb} \\
        -b $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.npa
    touch ${prefix}.npc
    touch ${prefix}.npl
    touch ${prefix}.npo

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' )
    END_VERSIONS
    """
}
