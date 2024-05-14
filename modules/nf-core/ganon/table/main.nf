process GANON_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::ganon=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ganon:1.5.1--py310h8abeb55_0':
        'biocontainers/ganon:1.5.1--py310h8abeb55_0' }"

    input:
    tuple val(meta), path(tre)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ganon \\
        table \\
        --input ${tre} \\
        --output-file ${prefix}.txt \\
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
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ganon: \$(echo \$(ganon --version 2>1) | sed 's/.*ganon //g')
    END_VERSIONS
    """
}
