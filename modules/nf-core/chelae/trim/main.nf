// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CHELAE_TRIM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f9f53b14969037fd2ef16d0c23b86ce41ce891d450e21d1402ff5217a2a8e5b/data'
        : 'community.wave.seqera.io/library/chelae:0.2.1--226b7d4244683014'}"

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
        --inputs ${reads} \\
        --outputs ${output_files} \\
        --threads ${task.cpus} \\
        --json ${prefix}.chelae.json \\
        --metrics ${prefix}.chelae.tsv \\
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
