process MACSYFINDER_DOWNLOAD {
    tag "${model_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macsyfinder:2.1.6--pyhdfd78af_0' :
        'biocontainers/macsyfinder:2.1.6--pyhdfd78af_0' }"

    input:
    val model_name

    output:
    path "models"       , emit: models
    tuple val("${task.process}"), val('macsyfinder'), eval('macsyfinder --version 2>&1 | sed "1!d;s/^.*MacSyFinder //;s/ .*$//"'), topic: versions, emit: versions_macsyfinder
    tuple val("${task.process}"), val('msf_data'), eval('msf_data --version 2>&1 | awk "NR==4" | sed "s/^- MacSyLib //;s/ *$//"'), topic: versions, emit: versions_macsydata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # msf_data installs models into the current directory by default
    # We'll create a models directory to store them
    mkdir -p models

    msf_data install \\
        --target models \\
        ${args} \\
        ${model_name}
    """

    stub:
    """
    mkdir -p models/${model_name}
    touch models/${model_name}/definitions.txt
    """
}
