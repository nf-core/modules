
process FAMSA_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::famsa=2.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/famsa:2.2.2--h9f5acd7_0':
        'biocontainers/famsa:2.2.2--h9f5acd7_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    famsa -gt_export \\
        $args \\
        -t ${task.cpus} \\
        ${fasta} \\
        ${prefix}.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$( famsa -help 2>&1 | head -n 2 | tail -n 1 | sed 's/ version //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$( famsa -help 2>&1 | head -n 2 | tail -n 1 | sed 's/ version //g' )
    END_VERSIONS
    """
}

