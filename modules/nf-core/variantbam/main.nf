process VARIANTBAM {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/variantbam:1.4.4a--h7d7f7ad_5' :
        'biocontainers/variantbam:1.4.4a--h7d7f7ad_5' }"

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
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.4.4a' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        variantbam: $VERSION
    END_VERSIONS
    """
}
