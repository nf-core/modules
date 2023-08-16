process GANON_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::ganon=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ganon:1.5.1--py310h8abeb55_0':
        'biocontainers/ganon:1.5.1--py310h8abeb55_0' }"

    input:
    tuple val(meta), path(rep)
    path db

    output:
    tuple val(meta), path("*.tre"), emit: tre
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    dbprefix=\$(find -L . -name '*.ibf' | sed 's/\\.ibf\$//')

    ganon \\
        report \\
        --input ${rep} \\
        --output-prefix ${prefix} \\
        --db-prefix \${dbprefix%%.ibf} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """
}
