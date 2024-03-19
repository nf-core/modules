process SEQFU_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.20.3--h1eb128b_0':
        'biocontainers/seqfu:1.20.3--h1eb128b_0' }"


    input:
    // stats can get one or more fasta or fastq files
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.tsv")    ,  emit: stats
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqfu \\
        stats \\
        $args \\
        --multiqc ${prefix}_mqc.txt \\
        $files > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mqc.txt
    seqfu stats ${prefix}_mqc.txt > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
