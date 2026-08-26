process HYPHY_GARD {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hyphy:2.5.101--h526e2cb_0':
        'quay.io/biocontainers/hyphy:2.5.101--h526e2cb_0' }"

    input:
    tuple val(meta), path(alignment)

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
    hyphy gard \\
        CPU=${task.cpus} \\
        --alignment ${alignment} \\
        --output ${prefix}_GARD.json \\
        ${args} \\
        > ${prefix}_GARD_output.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_GARD.json
    touch ${prefix}_GARD_output.txt
    """
}
