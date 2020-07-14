process INTERSECT_BED {
    tag "$input_file_1-$input_file_2"

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
    path (input_file_1)
    path (input_file_2)
    val (intersectbed_args)

    output:
    stdout()

    script:
    """
    bedtools intersect -a ${input_file_1} -b ${input_file_2} ${intersectbed_args}
    """
}

