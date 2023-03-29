process THERMORAWFILEPARSER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::thermorawfileparser=1.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/thermorawfileparser%3A1.4.2--ha8f3691_0' :
        'quay.io/biocontainers/thermorawfileparser:1.4.2--ha8f3691_0'    }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.mzML"), emit: mzML

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    thermorawfileparser \\
        -i=$raw \\
        -f=2 \\
        -o=./ \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echo \$(thermorawfileparser --version 2>&1)
    END_VERSIONS
    """
}
