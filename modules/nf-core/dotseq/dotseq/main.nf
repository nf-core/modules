process DOTSEQ_DOTSEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dotseq:1.0.0--r45hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-dotseq:1.0.0--r45hdfd78af_0' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts), path(flattened_gtf), path(flattened_bed)

    output:
    tuple val(meta), path("*.dou.interaction.dotseq.results.tsv"), emit: dou_interaction
    tuple val(meta), path("*.dte.interaction.dotseq.results.tsv"), emit: dte_interaction
    tuple val(meta), path("*.dou.strategy.dotseq.results.tsv")   , emit: dou_strategy   , optional: true
    tuple val(meta), path("*.dte.strategy.dotseq.results.tsv")   , emit: dte_strategy   , optional: true
    tuple val(meta), path("*.DOTSeqDataSets.rds")                , emit: rdata
    tuple val(meta), path("*.R_sessionInfo.log")                 , emit: session_info
    path "versions.yml"                                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dotseq.R'
}
