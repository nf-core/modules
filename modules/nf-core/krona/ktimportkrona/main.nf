process KRONA_KTIMPORTKRONA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    path html

    output:
    path "*.html"              , emit: html
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output = task.ext.prefix ? "${task.ext.prefix}.html" : 'krona.krona.html'
    """
    ktImportKrona ${html} \\
        -o ${output} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(ktImportKrona | grep -Po "(?<=KronaTools )[0-9.]+")
    END_VERSIONS
    """

    stub:
    def output = task.ext.prefix ? "${task.ext.prefix}.html" : 'krona.krona.html'
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$(ktImportKrona | grep -Po "(?<=KronaTools )[0-9.]+")
    END_VERSIONS
    """
}
