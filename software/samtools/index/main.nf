process samtools_index {
    tag "${bam.baseName}"

    container 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

    input:
    path(bam)

    output:
    path "*.bai"

    script:
    def suff_mem = ("${(task.memory.toBytes() - 6000000000) / task.cpus}" > 2000000000) ? 'true' : 'false'
    def avail_mem = (task.memory && suff_mem) ? "-m" + "${(task.memory.toBytes() - 6000000000) / task.cpus}" : ''
    """
    samtools index $bam \\
        -@ ${task.cpus} ${avail_mem}

    samtools --version &> v_samtools.txt
    """
}
