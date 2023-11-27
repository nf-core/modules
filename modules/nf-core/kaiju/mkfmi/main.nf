process KAIJU_MKFMI {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kaiju=1.9.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.9.2--h5b5514e_0':
        'biocontainers/kaiju:1.9.2--h5b5514e_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fmi"), emit: fmi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kaiju-mkbwt \\
        $args \\
        -n $task.cpus \\
        -o ${prefix} \\
        ${fasta}
    kaiju-mkfmi ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
