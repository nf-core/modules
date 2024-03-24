process BISCOT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscot:2.3.3--pyh7cba7a3_0':
        'biocontainers/biscot:2.3.3--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(cmap_ref), path(cmap1), path(cmap2), path(xmap1), path(xmap2), path(key), path(contigs), path(xmap_both)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.agp")  , emit: agp
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def xmap_both_arg = task.ext.args =~ "--xmap_2enz" ? "--xmap_2enz $xmap_both" : ""
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    biscot \\
        $args \\
        --cmap-ref $cmap_ref \\
        --cmap-1 $cmap1 \\
        --cmap-2 $cmap2 \\
        --xmap-1 $xmap1 \\
        --xmap-2 $xmap2 \\
        $xmap_both_arg \\
        --key $key \\
        --contigs $contigs \\
        --output .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscot: 2.3.3
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscot: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
