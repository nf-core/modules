process cutadapt {
    tag "${sample_id}"

    container 'quay.io/biocontainers/cutadapt:1.16--py27_1'

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple sample_id, file("trimmed_*.fastq")

    script:
    forward_fq = "trimmed_1.fastq"
    reverse_fq = "trimmed_2.fastq"


    if (params.singleEnd) {
        processing = """
                    cutadapt \
                        -j ${task.cpus} \
                        -q $params.cutadapt_min_quality \
                        --minimum-length $params.cutadapt_min_length \
                        --output ${forward_fq} \
                        ${reads}
                    """
    } else {
        processing = """
                    cutadapt \
                        -j ${task.cpus} \
                        -q $params.cutadapt_min_quality \
                        --minimum-length $params.cutadapt_min_length \
                        --pair-filter=any \
                        --output ${forward_fq} \
                        --paired-output ${reverse_fq} ${reads}


                    """
    }

    version = """
    cutadapt --version &> v_cutadapt.txt
    """

    return processing + version
}
