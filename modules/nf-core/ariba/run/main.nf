process ARIBA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_3':
        'biocontainers/ariba:2.14.6--py39h67e14b5_3' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    tuple val("${task.process}"), val('ariba'), eval('ariba version 2>/dev/null | head -1 | sed "s/ARIBA version: //"'), emit: versions_ariba, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    ariba \\
        run \\
        ${db_name}/ \\
        ${reads} \\
        ${prefix} \\
        ${args} \\
        --threads ${task.cpus}
    """
}
