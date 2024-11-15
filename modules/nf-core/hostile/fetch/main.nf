process HOSTILE_FETCH {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hostile:1.1.0--pyhdfd78af_0':
        'biocontainers/hostile:1.1.0--pyhdfd78af_0' }"

    output:
    path "reference/"   , emit: reference
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir reference/
    export HOSTILE_CACHE_DIR=./reference

    hostile \\
        fetch \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir reference/
    export HOSTILE_CACHE_DIR=./reference

    touch reference/human-t2t-hla.1.bt2
    touch reference/human-t2t-hla.2.bt2
    touch reference/human-t2t-hla.3.bt2
    touch reference/human-t2t-hla.4.bt2
    touch reference/human-t2t-hla.rev.1.bt2
    touch reference/human-t2t-hla.rev.2.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """
}
