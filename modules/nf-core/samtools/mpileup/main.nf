process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0' :
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(input)       // sam/bam/cram
    tuple val(meta), path(index)       // csi/bai/crai - optional, only required when using '-r'
    tuple val(meta), path(intervals)   // optional
    tuple val(meta2), path(fasta)      // optional

    output:
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed -e "s/samtools //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def fasta_cmd     = fasta           ? "--fasta-ref $fasta" : ""
    def intervals_cmd = intervals       ? "-l ${intervals}" : ""
    """
    samtools mpileup \\
        $fasta_cmd \\
        --output ${prefix}.mpileup \\
        $args \\
        $intervals_cmd \\
        $input
    bgzip ${prefix}.mpileup
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.mpileup.gz
    """
}
