process WISECONDORX_GENDER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), stdout , emit: gender
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    WisecondorX gender \\
        ${npz} \\
        ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def VERSION = '1.2.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo male

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
