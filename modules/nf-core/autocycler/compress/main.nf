
process AUTOCYCLER_COMPRESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/autocycler:0.5.2--h3ab6199_0':
        'biocontainers/autocycler:0.5.2--h3ab6199_0' }"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("$prefix/*.gfa"),  emit: gfa
    tuple val(meta), path("$prefix/*.yaml"), emit: stats
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    autocycler compress \\
        $args \\
        -t $task.cpus \\
        -i . \\
        -a ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    touch ${prefix}/input_assemblies.gfa
    touch ${prefix}/input_assemblies.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version |  sed 's/^[^ ]* //')
    END_VERSIONS
    """
}
