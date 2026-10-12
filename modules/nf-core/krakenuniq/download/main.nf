process KRAKENUNIQ_DOWNLOAD {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2'
        : 'quay.io/biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2'}"

    input:
    val pattern

    output:
    path "${pattern}/", emit: output
    tuple val("${task.process}"), val('krakenuniq'), eval("krakenuniq --version | sed '1!d;s/KrakenUniq version //'"), emit: versions_krakenuniq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    krakenuniq-download \\
        --threads ${task.cpus} \\
        -o ${pattern}/ \\
        ${pattern} \\
        ${args}
    """

    stub:
    """
    mkdir ${pattern}
    """
}
