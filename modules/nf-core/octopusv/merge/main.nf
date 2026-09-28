process OCTOPUSV_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8378c503810508764badf4504647ea1791fa1f585dc5b77ad5710237609b15ff/data'
:         'community.wave.seqera.io/library/octopusv:1.0.0--5a8c89c65098e801' }"

    input:
    tuple val(meta), path(svcfs, arity: '1..*'), path(input_list), path(specific_svcfs, arity: '0..*'), val(strategy_flag)

    output:
    tuple val(meta), path("*.svcf"), emit: svcf
    tuple val(meta), path("*.png"), emit: upset_plot, optional: true
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

    def input_files_arg = input_list ? '' : svcfs.collect { svcf -> "--input-file ${svcf}" }.join(' ')
    def input_list_arg = input_list ? "--input-list ${input_list}" : ''
    def merge_strategy_arg = merge_strategy == 'specific'
        ? specific_svcfs.collect { svcf -> "--specific ${svcf}" }.join(' ')
        : "--${merge_strategy}"

    """
    octopusv merge ${input_files_arg} \\
        ${input_list_arg} \\
        --output-file ${prefix}.svcf \\
        ${merge_strategy_arg} \\
        ${args}
    """

    stub:
    def args = (task.ext.args ?: '').trim()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_strategy = (strategy_flag ?: '').trim()
    if (!merge_strategy) {
        merge_strategy = 'union'
    }

    def input_files_arg = input_list ? '' : svcfs.collect { svcf -> "--input-file ${svcf}" }.join(' ')
    def input_list_arg = input_list ? "--input-list ${input_list}" : ''
    def merge_strategy_arg = merge_strategy == 'specific'
        ? specific_svcfs.collect { svcf -> "--specific ${svcf}" }.join(' ')
        : "--${merge_strategy}"
    """
    echo ${input_files_arg} ${input_list_arg} ${merge_strategy_arg} ${args}

    touch ${prefix}.svcf
    """
}
