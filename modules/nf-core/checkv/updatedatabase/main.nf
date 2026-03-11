process CHECKV_UPDATEDATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0':
        'biocontainers/checkv:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path (fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/*"), emit: checkv_db
    tuple val("${task.process}"), val("checkv"), eval("checkv -h 2>&1 | sed '1!d;s/^.*CheckV v//;s/:.*//'"), topic: versions, emit: versions_checkv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def checkv_db = db ?: ''
    def update_sequence = fasta ?: ''
    """
    checkv update_database \\
        --threads $task.cpus \\
        $args \\
        $checkv_db \\
        ./$prefix/  \\
        $update_sequence
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/
    touch ${prefix}/README.txt
    """

}
