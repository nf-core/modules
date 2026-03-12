process CLUSTALO_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5':
        'biocontainers/clustalo:1.2.4--h87f3376_5' }"

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
    clustalo \\
        -i ${fasta} \\
        --guidetree-out ${prefix}.dnd \\
        --threads=${task.cpus} \\
        $args

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd

    cat <<-END_VERSIONS > tuple val("${task.process}"), val('<tool1>'), eval('tool1 --version'), emit: versions_tool1, topic: versions
    "${task.process}":
        clustalo: \$( clustalo --version )
    END_VERSIONS
    """
}
