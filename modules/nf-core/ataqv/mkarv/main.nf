process ATAQV_MKARV {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ataqv:1.3.1--py310ha155cf9_1':
        'biocontainers/ataqv:1.3.1--py310ha155cf9_1' }"

    input:
    path "jsons/*"

    output:
    path "html"        , emit: html
    tuple val("${task.process}"), val('ataqv'), eval('echo \$(ataqv --version)'), emit: versions_ataqv, topic: versions
    // tuple val("${task.process}"), val('mkarv'), eval('mkarv --version'), emit: versions_mkarv, topic: versions //Use this when version string has been fixed
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkarv \\
        ${args} \\
        --concurrency ${task.cpus} \\
        --force \\
        ./html/ \\
        jsons/*
    """

    stub:
    """
    mkdir -p html
    touch html/index.html
    """
}
