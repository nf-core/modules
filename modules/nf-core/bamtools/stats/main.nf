process BAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bamtools=2.5.1" : null)
    def container_image = "/bamtools:2.5.1--h9a82719_9"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bamtools \\
        stats \\
        -in $bam \\
        >${prefix}.bam.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
