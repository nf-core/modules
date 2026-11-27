process REFINEGEMS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.16' :
        'docker.io/library/python:3.10-slim' }"

    input:
    tuple val(meta), val(config_type)

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    refinegems setup config \\
        --type ${config_type} \\
        --filename ${prefix}.yaml \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        refinegems: \$(refinegems --version | sed -E 's/.*([0-9]+\.[0-9]+\.[0-9]+(\.[0-9]+)?).*/\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        refinegems: 0.1.23
    END_VERSIONS
    """
}
