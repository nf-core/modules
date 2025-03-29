process DREP_COMPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0'
        : 'biocontainers/drep:3.5.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${prefix}"), emit: directory
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    dRep \\
        compare \\
        ${prefix} \\
        -p ${task.cpus} \\
        ${args} \\
        -g fastas/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep | head -n 2 | sed 's/.*v//g;s/ .*//g' | tail -n 1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep | head -n 2 | sed 's/.*v//g;s/ .*//g' | tail -n 1)
    END_VERSIONS
    """
}
