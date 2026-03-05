process BIGSLICE {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided correctly by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0':
        'biocontainers/bigslice:2.0.2--pyh8ed023e_0' }"

    input:
    tuple val(meta), path(bgc)
    path(hmmdb)

    output:
    tuple val(meta), path("${prefix}/result/data.db")    , emit: db
    tuple val(meta), path("${prefix}/result/tmp/**/*.fa"), emit: fa
    tuple val("${task.process}"), val('bigslice'), eval("echo 2.0.2"), topic: versions, emit: versions_bigslice

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bigslice \\
        $args \\
        --num_threads ${task.cpus} \\
        -i ${bgc} \\
        --program_db_folder ${hmmdb} \\
        ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}/result/tmp/2e555308dfc411186cf012334262f127
    touch ${prefix}/result/data.db
    touch ${prefix}/result/tmp/2e555308dfc411186cf012334262f127/test.fa
    """
}
