process HOMER_MAKETAGDIRECTORY {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0fe4a3875b78dce3c66b43fb96489769cc32e55e329e2525d2af09096af2252a/data'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-edger_homer_samtools_pruned:a8f4c58755bb281b'}"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_tagdir"), emit: tagdir
    tuple val(meta), path("*_tagdir/tagInfo.txt"), emit: taginfo
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('homer'), val("4.11"), emit: versions_homer, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | sed '1!d;s/^.*samtools //'"), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('deseq2'), eval('Rscript -e "cat(as.character(packageVersion(\'DESeq2\')))"'), emit: versions_deseq2, topic: versions
    tuple val("${task.process}"), val('edger'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'), emit: versions_edger, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeTagDirectory \\
        ${prefix}_tagdir \\
        -genome ${fasta} \\
        ${args} \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}_tagdir/

    touch ${prefix}_tagdir/${fasta.baseName}.1.tags.tsv
    touch ${prefix}_tagdir/tagAutocorrelation.txt
    touch ${prefix}_tagdir/tagCountDistribution.txt
    touch ${prefix}_tagdir/tagInfo.txt
    touch ${prefix}_tagdir/tagLengthDistribution.txt
    """
}
