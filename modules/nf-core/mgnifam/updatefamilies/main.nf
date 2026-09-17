process MGNIFAM_UPDATEFAMILIES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnifam:3.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/mgnifam:3.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(hmms, stageAs: 'hmm_input/*'), path(fasta_file), path(fasta_index)

    output:
    tuple val(meta), path("${prefix}/seed_msa/*.sto.gz")                , emit: seed_msa  , optional: true
    tuple val(meta), path("${prefix}/full_msa/*.sto.gz")                , emit: full_msa  , optional: true
    tuple val(meta), path("${prefix}/hmm/*.hmm.gz")                     , emit: hmm       , optional: true
    tuple val(meta), path("${prefix}/rf/*.txt")                         , emit: rf        , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_families.tsv")   , emit: tsv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_metadata.csv")   , emit: csv       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated.log")            , emit: log       , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_reps.fasta.gz")  , emit: reps_fasta, optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_successful.txt") , emit: successful, optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_discarded.csv")  , emit: discarded , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_converged.txt")  , emit: converged , optional: true
    tuple val(meta), path("${prefix}/${prefix}_updated_delta.csv")      , emit: delta     , optional: true
    tuple val("${task.process}"), val('mgnifam'), eval("mgnifam --version 2>&1"), topic: versions, emit: versions_mgnifam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def index = fasta_index ? "--fasta_index ${fasta_index}" : ''
    """
    lib_file=\$(find hmm_input -maxdepth 1 -name '*.hmm.lib*' -print -quit)

    mgnifam update_families \\
        --hmm_input "\${lib_file:-hmm_input}" \\
        --fasta_file "${fasta_file}" \
        --output_dir "${prefix}" \
        --cpus ${task.cpus} \
        --chunk_id "${prefix}" \
        ${index} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir -p ${prefix}/seed_msa ${prefix}/full_msa ${prefix}/hmm ${prefix}/rf
    echo "" | gzip > ${prefix}/seed_msa/1_7.sto.gz
    echo "" | gzip > ${prefix}/full_msa/1_7.sto.gz
    echo "" | gzip > ${prefix}/hmm/1_7.hmm.gz
    touch ${prefix}/rf/1_7.txt
    touch ${prefix}/${prefix}_updated_families.tsv
    touch ${prefix}/${prefix}_updated_metadata.csv
    touch ${prefix}/${prefix}_updated.log
    echo "" | gzip > ${prefix}/${prefix}_updated_reps.fasta.gz
    touch ${prefix}/${prefix}_updated_successful.txt
    touch ${prefix}/${prefix}_updated_discarded.csv
    touch ${prefix}/${prefix}_updated_converged.txt
    touch ${prefix}/${prefix}_updated_delta.csv
    """
}
