process ARIBA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ariba=2.14.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_3':
        'biocontainers/ariba:2.14.6--py39h67e14b5_3' }"

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    ariba \\
        run \\
        ${db_name}/ \\
        ${reads} \\
        ${prefix} \\
        $args \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
