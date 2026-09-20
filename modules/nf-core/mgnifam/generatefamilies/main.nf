process MGNIFAM_GENERATEFAMILIES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnifam:3.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/mgnifam:3.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(clustering), path(fasta_file), path(fasta_index)

    output:
    tuple val(meta), path("${prefix}/seed_msa/*.sto.gz")       , emit: seed_msa  , optional: true
    tuple val(meta), path("${prefix}/full_msa/*.sto.gz")       , emit: full_msa  , optional: true
    tuple val(meta), path("${prefix}/hmm/*.hmm.gz")            , emit: hmm       , optional: true
    tuple val(meta), path("${prefix}/rf/*.txt")                , emit: rf        , optional: true
    tuple val(meta), path("${prefix}/${prefix}_families.tsv")  , emit: tsv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_metadata.csv")  , emit: csv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}.log")           , emit: log       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_reps.fasta.gz") , emit: reps_fasta, optional: true
    tuple val(meta), path("${prefix}/${prefix}_successful.txt"), emit: successful, optional: true
    tuple val(meta), path("${prefix}/${prefix}_discarded.csv") , emit: discarded , optional: true
    tuple val(meta), path("${prefix}/${prefix}_converged.txt") , emit: converged , optional: true
    tuple val("${task.process}"), val('mgnifam'), eval("mgnifam --version 2>&1"), topic: versions, emit: versions_mgnifam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def index = fasta_index ? "--fasta_index ${fasta_index}" : ''
    """
    mgnifam generate_families \\
        --clusters_chunk ${clustering} \\
        --fasta_file ${fasta_file} \\
        --output_dir ${prefix} \\
        --cpus ${task.cpus} \\
        --chunk_id ${prefix} \\
        ${index} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}/seed_msa ${prefix}/full_msa ${prefix}/hmm ${prefix}/rf
    echo "" | gzip > ${prefix}/seed_msa/${prefix}_1.sto.gz
    echo "" | gzip > ${prefix}/full_msa/${prefix}_1.sto.gz
    echo "" | gzip > ${prefix}/hmm/${prefix}_1.hmm.gz
    touch ${prefix}/rf/${prefix}_1.txt
    touch ${prefix}/${prefix}_families.tsv
    touch ${prefix}/${prefix}_metadata.csv
    touch ${prefix}/${prefix}.log
    echo "" | gzip > ${prefix}/${prefix}_reps.fasta.gz
    touch ${prefix}/${prefix}_successful.txt
    touch ${prefix}/${prefix}_discarded.csv
    touch ${prefix}/${prefix}_converged.txt
    """
}
