process gatk_dict {
    tag {fasta}

    container 'quay.io/biocontainers/gatk4-spark:4.1.4.1--1'

    input:
        path fasta

    output:
        path "${fasta.baseName}.dict"

    script:
    """
    gatk --java-options "-Xmx${task.memory.giga}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}
