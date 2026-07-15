process CHELAE_TRIM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/chelae:0.1.0--hfa8f182_0'
        : 'quay.io/biocontainers/chelae:0.1.0--hfa8f182_0'}"

    input:
    tuple val(meta), path(reads)
    path adapter_fasta

    output:
    tuple val(meta), path('*.chelae.fastq.gz'), emit: reads
    tuple val(meta), path('*.chelae.json'), emit: json
    tuple val(meta), path('*.chelae.tsv'), emit: metrics
    tuple val("${task.process}"), val('chelae'), eval("chelae --version 2>&1 | sed 's/chelae //'"), emit: versions_chelae, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_fasta_arg = adapter_fasta ? "--adapter-fasta ${adapter_fasta}" : ''
    def output_files = meta.single_end ? "${prefix}.chelae.fastq.gz" : "${prefix}_R1.chelae.fastq.gz ${prefix}_R2.chelae.fastq.gz"
    """
    chelae trim \\
        -i ${reads} \\
        -o ${output_files} \\
        -t ${task.cpus} \\
        -j ${prefix}.chelae.json \\
        -m ${prefix}.chelae.tsv \\
        ${adapter_fasta_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_reads = meta.single_end
        ? "echo '' | gzip > ${prefix}.chelae.fastq.gz"
        : "echo '' | gzip > ${prefix}_R1.chelae.fastq.gz ; echo '' | gzip > ${prefix}_R2.chelae.fastq.gz"
    """
    ${create_reads}
    echo '{}' > ${prefix}.chelae.json
    printf 'metric\tvalue\n' > ${prefix}.chelae.tsv
    """
}
