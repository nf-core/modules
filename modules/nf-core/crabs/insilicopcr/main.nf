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
    tuple val(meta), path("${prefix}.txt"), emit: txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.insilicopcr"
    def version_cmd = "\$(crabs --help 2>/dev/null | grep 'CRABS |' | sed 's/.*CRABS | v\\([0-9.]*\\).*/\\1/')"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        version_cmd    = '1.0.7'
    }
    """
    crabs --in-silico-pcr \\
        --input ${crabsdb} \\
        --output ${prefix}.txt \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: ${version_cmd}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.insilicopcr"
    def version_cmd = "\$(crabs --help 2>/dev/null | grep 'CRABS |' | sed 's/.*CRABS | v\\([0-9.]*\\).*/\\1/')"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        version_cmd = '1.0.7'
    }
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crabs: ${version_cmd}
    END_VERSIONS
    """
}
