process BEDTOOLS_MERGE {
    tag { input_file }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
    path (input_file)
    val (bedtools_merge_args)

    output:
    stdout()

    script:
    """
    bedtools merge -i ${input_file} ${bedtools_merge_args}
    """
}
