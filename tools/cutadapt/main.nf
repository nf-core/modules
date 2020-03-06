process cutadapt {
    tag "${sample_id}"

    container 'quay.io/biocontainers/cutadapt:1.16--py27_1'

    input:
    tuple sample_id, file(input_forward_fq), file(input_reverse_fq)

    output:
    tuple sample_id, file(forward_fq), file(reverse_fq)

    script:
    forward_fq = "trimmed_forward.fastq"
    reverse_fq = "trimmed_reverse.fastq"

    """
    cutadapt \
			-j ${task.cpus} \
			-q $params.cutadapt_min_quality \
			--minimum-length $params.cutadapt_min_length \
			--pair-filter=any \
			--output ${forward_fq} \
			--paired-output ${reverse_fq} '${input_forward_fq}' '${input_reverse_fq}'
    
    cutadapt --version &> v_cutadapt.txt
    """
}
