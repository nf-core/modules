process RIBOTISH_QUALITY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotish:0.2.7--pyhdfd78af_0':
        'biocontainers/ribotish:0.2.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt")    , emit: distribution
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.para.py"), emit: offset
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ribotish quality \\
        -b $bam \\
        -g $gtf \\
        -o ${prefix}_qual.txt \\
        -f ${prefix}_qual.pdf \\
        -r ${prefix}.para.py \\
        -p $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """

//    stub:
//    def args = task.ext.args ?: ''
//    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
//    """
//    touch ${prefix}.bam

//    cat <<-END_VERSIONS > versions.yml
//    "${task.process}":
//        ribotish: \$(samtools --version |& sed '1!d ; s/samtools //')
//    END_VERSIONS
//    """
}
