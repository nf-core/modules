process WISECONDORX_GENDER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::wisecondorx=1.2.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.5--pyh5e36f6f_0':
        'quay.io/biocontainers/wisecondorx:1.2.5--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(npz)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), stdout , emit: gender
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

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
    def VERSION = '1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    echo male

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wisecondorx: ${VERSION}
    END_VERSIONS
    """
}
