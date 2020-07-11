process samtools_faidx {
    tag {fasta}

    container 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

    input:
        path(fasta)

    output:
        path("${fasta}.fai")

    script:
    """
    samtools faidx ${fasta}
    """
}
