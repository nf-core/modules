process CLEANIFIER_DOWNLOAD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cleanifier:1.3.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/cleanifier:1.3.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(version)

    output:
    tuple val(meta), path("*.{filter,hash}"), emit: index
    tuple val(meta), path("*.info"), emit: info
    tuple val("${task.process}"), val('cleanifier'), eval("cleanifier --version"), emit: versions_cleanifier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cleanifier download \\
        --version ${version} \\
        --dir . \\
        ${args}

    rm *.tar
    """

    stub:
    def suffix = (version == "exact") ? "hash" : "filter"
    """
    touch index.${suffix}
    touch index.info
    """
}
