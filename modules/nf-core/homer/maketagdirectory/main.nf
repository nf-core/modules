
process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0' :
        'biocontainers/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0' }"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_tagdir")            , emit: tagdir
    tuple val(meta), path("*_tagdir/tagInfo.txt"), emit: taginfo
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    makeTagDirectory \\
        ${prefix}_tagdir \\
        -genome $fasta \\
        $args \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
