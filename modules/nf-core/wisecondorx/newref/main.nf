process WISECONDORX_NEWREF {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::wisecondorx=1.2.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.5--pyh5e36f6f_0':
        'biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    inputs.each { if("${it}" == "${prefix}.npz") error "${it} has the same name as the output file, set prefix in module configuration to disambiguate!"}

    """
    WisecondorX \\
        newref \\
        *.npz \\
        ${prefix}.npz \\
        --cpus ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    inputs.each { if("${it}" == "${prefix}.npz") error "${it} has the same name as the output file, set prefix in module configuration to disambiguate!"}

    """
    touch ${prefix}.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
