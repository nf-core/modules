process VARIANTBAM {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::variantbam=1.4.4a" : null)
        def container_image = "/variantbam:1.4.4a--h7d7f7ad_5"
                                             container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.4.4a' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    variant \\
        $bam \\
        -o ${prefix}.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        variantbam: $VERSION
    END_VERSIONS
    """
}
