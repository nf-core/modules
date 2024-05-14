process CLUSTALO_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::clustalo=1.2.4"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustalo: \$( clustalo --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustalo: \$( clustalo --version )
    END_VERSIONS
    """
}
