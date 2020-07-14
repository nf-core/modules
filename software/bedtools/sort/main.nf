process BEDTOOLS_SORT {
    tag { input_file }

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
    path (input_file)
    val (bedtools_sort_args)

    output:
    stdout()

    script:
    """
    bedtools sort -i ${input_file} ${bedtools_sort_args}
    """
}
