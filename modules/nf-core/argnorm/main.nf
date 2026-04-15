process ARGNORM {
    tag "${meta.id}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/argnorm:0.8.0--pyhdfd78af_0'
        : 'biocontainers/argnorm:0.8.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(input_tsv)
    val tool
    val db

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('argnorm'), eval('argnorm --version'), emit: versions_argnorm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def db_args = db ? "--db ${db}" : ""
    if (!tool) {
        error('Tool not provided.')
    }
    if ((tool in ["abricate"]) && !db) {
        error("${tool} requires a database but <db> not provided.")
    }

    """
    argnorm \\
        ${tool} \\
        -i ${input_tsv} \\
        -o ${prefix} \\
        ${db_args} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!tool) {
        error('Tool not provided.')
    }
    if ((tool in ["abricate"]) && !db) {
        error("${tool} requires a database but <db> not provided.")
    }

    """
    touch ${prefix}
    """
}
