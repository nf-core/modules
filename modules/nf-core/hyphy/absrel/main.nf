process HYPHY_ABSREL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hyphy:2.5.101--h526e2cb_0':
        'quay.io/biocontainers/hyphy:2.5.101--h526e2cb_0' }"

    input:
    tuple val(meta), path(alignment), path(tree)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*_output.txt"), emit: log
    tuple val("${task.process}"), val('hyphy'), eval("hyphy --version | sed -E 's/.*HYPHY ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/'"), topic: versions, emit: versions_hyphy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hyphy absrel \\
        CPU=${task.cpus} \\
        --alignment ${alignment} \\
        --tree ${tree} \\
        --output ${prefix}_ABSREL.json \\
        ${args} \\
        > ${prefix}_ABSREL_output.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_ABSREL.json
    touch ${prefix}_ABSREL_output.txt
    """
}
