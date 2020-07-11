process bwa_index {
    tag {fasta}

    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    input:
        path(fasta)

    output:
        path("${fasta}.*")

    script:
    """
    bwa index ${fasta}
    """
}