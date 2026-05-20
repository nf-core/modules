process OCTOPUSV_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/octopusv:0.3.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/octopusv:0.3.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(svcfs), val(strategy_flag)

    output:
    tuple val(meta), path("*.svcf"), emit: svcf
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = (task.ext.args ?: '').trim()
    def prefix = task.ext.prefix ?: "${meta.id}"

    def merge_strategy = (strategy_flag ?: '').trim()
    if (!merge_strategy) {
        merge_strategy = 'union'
    }

    """
    octopusv merge -i ${svcfs} \\
        -o ${prefix}.svcf \\
        --${merge_strategy} \\
        ${args}
    """

    stub:
    def args = (task.ext.args ?: '').trim()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_strategy = (strategy_flag ?: '').trim()
    if (!merge_strategy) {
        merge_strategy = '--union'
    }
    """
    echo ${merge_strategy} ${args}

    touch ${prefix}.svcf
    """
}
