process NONPAREIL_NONPAREIL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.4.1--r42h4ac6f70_4':
        'biocontainers/nonpareil:3.4.1--r42h4ac6f70_4' }"

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
    def is_compressed = reads.getExtension() == "gz" ? true : false
    def reads_name = is_compressed ? reads.getBaseName() : reads
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${reads} > ${reads_name}
    fi

    nonpareil \\
        -s $reads_name \\
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
