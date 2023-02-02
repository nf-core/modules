process SAMBAMBA_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sambamba=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'maulik23/sambamba:1.0':
        'maulik23/sambamba:1.0' }"

    input:
        tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.mrkDup.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sambamba \\
        markdup \\
        $args \\
        -t $task.cpus \\
        $bam \\
        ${prefix}.mrkDup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}
