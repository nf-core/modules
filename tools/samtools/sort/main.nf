process samtools_index {
    tag "${bam.baseName}"

    container: 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'

    input:
    set file(bam)

    output:
    file "*.bam.bai"

    script:
    """
    samtools index $bam

    samtools --version &> v_samtools.txt
    """
}
