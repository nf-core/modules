process BIGSLICE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0':
        'biocontainers/bigslice:2.0.2--pyh8ed023e_0' }"

    input:
    tuple val(meta), path(bgc)
    path(hmmdb)

    output:
    tuple val(meta), path("${prefix}/result/data.db")    , emit: db
    // tuple val(meta), path("${prefix}/result/tmp/**/*.fa"), emit: fa
    tuple val("${task.process}"), val('bigslice'), eval("bigslice --version"), topic: versions, emit: versions_bigslice

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bigslice \\
        $args \\
        -i ${bgc} \\
        --program_db_folder ${hmmdb} \\
        ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.bam
    """
}
