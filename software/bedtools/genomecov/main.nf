process BEDTOOLS_GENOMECOV {
    tag {bam}

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
    path (bam)
    path (chrom_sizes)
    val (bedtools_genomecov_args)

    output:
    stdout()

    script:
    """
    bedtools genomecov -ibam ${bam} -g ${chrom_sizes} ${bedtools_genomecov_args}
    """
}
