process CRABS_INSILICOPCR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crabs:1.0.7--pyhdfd78af_0':
        'biocontainers/crabs:1.0.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(crabsdb)

    output:
    tuple val(meta), path("*.insilicopc.txt"), emit: txt
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    crabs --in-silico-pcr \\
        --input ${crabsdb} \\
        --output ${prefix}.insilicopcr.txt \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --help | grep 'CRABS |' | sed 's/.*CRABS | \\(v[0-9.]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.insilicopcr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: \$(crabs --help | grep 'CRABS |' | sed 's/.*CRABS | \\(v[0-9.]*\\).*/\\1/')
    END_VERSIONS
    """
}
