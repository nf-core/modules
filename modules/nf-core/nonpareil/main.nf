
process NONPAREIL {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::nonpareil=3.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.4.1--r42h9f5acd7_3':
        'biocontainers/nonpareil:3.4.1--r42h9f5acd7_3' }"

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
        : \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' ))
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
        : \$(echo \$(nonpareil -V 2>&1) | sed 's/Nonpareil v//' ))
    END_VERSIONS
    """
}
