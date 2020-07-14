process GENOMECOV {
    tag {bam}

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        path (bam)
        path (fasta_sizes)
        val (genomecov_args)

    output:
        stdout()

    script:
    """
    bedtools genomecov -ibam ${bam} -g ${fasta_sizes} ${genomecov_args}
    """
}
