process HOMER_MAKETAGDIRECTORY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // singularity build url: https://wave.seqera.io/view/builds/bd-9c603739ae7d4fd3_1
    // docker build url: https://wave.seqera.io/view/builds/bd-08c7bb832e96c6bd_1
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:9c603739ae7d4fd3'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:08c7bb832e96c6bd'}"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_tagdir"), emit: tagdir
    tuple val(meta), path("*_tagdir/tagInfo.txt"), emit: taginfo
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    makeTagDirectory \\
        ${prefix}_tagdir \\
        -genome ${fasta} \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '5.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir ${prefix}_tagdir/

    touch ${prefix}_tagdir/${fasta.baseName}.1.tags.tsv
    touch ${prefix}_tagdir/tagAutocorrelation.txt
    touch ${prefix}_tagdir/tagCountDistribution.txt
    touch ${prefix}_tagdir/tagInfo.txt
    touch ${prefix}_tagdir/tagLengthDistribution.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
