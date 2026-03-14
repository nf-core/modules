process KRONA_KTIMPORTKRONA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1'
        : 'biocontainers/krona:2.8.1--pl5321hdfd78af_1'}"

    input:
    path html

    output:
    path "${prefix}.html", emit: html
    tuple val("${task.process}"), val('krona'), eval("ktImportKrona | grep -Po '(?<=KronaTools )[0-9.]+'"), topic: versions, emit: versions_krona

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ? "${task.ext.prefix}" : 'krona.krona'
    """
    ktImportKrona ${html} \\
        -o ${prefix}.html \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ? "${task.ext.prefix}" : 'krona.krona'
    """
    touch ${prefix}.html
    """
}
