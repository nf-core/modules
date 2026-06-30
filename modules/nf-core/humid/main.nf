
process HUMID {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humid:1.0.4--hadf994f_0':
        'quay.io/biocontainers/humid:1.0.4--hadf994f_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(umi_file)

    output:
    tuple val(meta), path("${prefix}.log")         , emit: log
    tuple val(meta), path("*_dedup*.fastq.gz")     , emit: dedup    , optional: true
    tuple val(meta), path("*_annotated*.fastq.gz") , emit: annotated, optional: true
    tuple val(meta), path("${prefix}")             , emit: stats    , optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('humid'), val("1.0.4"), emit: versions_humid, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def umis = umi_file ?: ''

    """
    humid \\
        ${args} \\
        -d ${prefix} \\
        -l ${prefix}.log \\
        ${reads} \\
        ${umis} \\

    mv ${prefix}/*.fastq* .
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo "" | gzip > ${prefix}_1_dedup.fastq.gz
    echo "" | gzip > ${prefix}_2_dedup.fastq.gz
    echo "" | gzip > ${prefix}_1_annotated.fastq.gz
    echo "" | gzip > ${prefix}_2_annotated.fastq.gz
    touch ${prefix}/stats.dat
    touch ${prefix}/neigh.dat
    touch ${prefix}/counts.dat
    touch ${prefix}/clusters.dat
    touch ${prefix}.log
    """
}
