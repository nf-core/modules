process CLUSTALO_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clustalo:1.2.4--h87f3376_5':
        'biocontainers/clustalo:1.2.4--h87f3376_5' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(tree)

    output:
    tuple val(meta), path("*.aln.gz"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clustalo \\
        -i ${fasta} \\
        --threads=${task.cpus} \\
        $args \\
        | pigz -p ${task.cpus} -c > ${prefix}.aln.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustalo: \$( clustalo --version )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustalo: \$( clustalo --version )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
