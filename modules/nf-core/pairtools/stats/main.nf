process PAIRTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0' :
        'biocontainers/pairtools:1.0.2--py39h2a9f597_0' }"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path("*.pairs.stat"), emit:stats
    path("versions.yml"), emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools stats \\
        ${args} \\
        --nproc-in ${task.cpus} --nproc-out ${task.cpus} \\
        -o ${prefix}.pairs.stat \\
        ${pairs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
}
