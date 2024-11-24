
process HOMER_MAKETAGDIRECTORY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/898620b4f7461c66fc41f8532d3ae1cb6ae9f44d304a59c988893244dbcbe389/data' :
        'community.wave.seqera.io/library/homer_wget_samtools_r-essentials_pruned:60f842302c4fb688' }"

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

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir ${prefix}_tagdir/

    touch ${prefix}_tagdir/${fasta.baseName}.1.tags.tsv
    touch ${prefix}_tagdir/tagAutocorrelation.txt
    touch ${prefix}_tagdir/tagCountDistribution.txt
    touch ${prefix}_tagdir/tagInfo.txt
    touch ${prefix}_tagdir/tagLengthDistribution.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
