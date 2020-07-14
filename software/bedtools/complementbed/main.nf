process COMPLEMENT_BED {
    tag {input_file}

    container 'quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0'

    input:
        path (input_file)
        path (fasta_sizes)
        val (complementbed_args)

    output:
        stdout()

    script:
    """
    bedtools complement -i ${input_file} -g ${fasta_sizes} ${complementbed_args}
    """
}
