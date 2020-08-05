process bwa-mem2_index {
    tag {fasta}

    container 'quay.io/biocontainers/bwa-mem2:2.0--he513fc3_0'

    input:
        path(fasta)

    output:
        path("${fasta}.*")

    script:
    """
    bwa-mem2 index ${fasta}
    """
}
