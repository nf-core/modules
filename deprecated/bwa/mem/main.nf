params.bwa_options = "-M -B 2"
params.sequencer = "ILLUMINA"

process bwa_mem {
    tag "$id"

    publishDir "${params.outdir}/bwa_mem", mode: 'copy'

    //TO-DO: Change container declaration, for now a test container is present in my personal docker acccount
    container 'jeremy1805/bwa-mem-img'

    input:
        tuple val(id), path(reads)
        path genomeindex
        val indexprefix

    output:
        tuple path("*.bam"), path("*.bai")

    script:
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${id}\\tLB:${id}\\tSM:${id}\\tPL:${params.sequencer}" \\
    ${params.bwa_options} ${indexprefix} ${reads} | samtools sort -@8 -O BAM -o ${id}.bam -

    samtools index ${id}.bam
    """
}
