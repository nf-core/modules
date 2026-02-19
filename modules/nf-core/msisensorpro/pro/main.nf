process MSISENSORPRO_PRO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msisensor-pro%3A1.3.0--hfef96ef_0':
        'biocontainers/msisensor-pro:1.3.0--hfef96ef_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(list)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("${prefix}")          , emit: summary_msi
    tuple val(meta), path("${prefix}_all")      , emit: all_msi
    tuple val(meta), path("${prefix}_dis")      , emit: dis_msi
    tuple val(meta), path("${prefix}_unstable") , emit: unstable_msi
    tuple val("${task.process}"), val('msisensor-pro'), eval("msisensor-pro --version 2>&1 | sed -nE 's/Version:\\s*//p'") , emit: versions_msisensorpro, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor-pro \\
        pro \\
        -d $list \\
        -t $input \\
        -g $fasta \\
        -o ${prefix} \\
        -b $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix} ${prefix}_all ${prefix}_dis ${prefix}_unstable
    """
}
