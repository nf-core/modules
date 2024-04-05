process ANOTA2SEQ_ANOTA2SEQRUN {
    tag "$meta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-anota2seq:1.24.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-anota2seq:1.24.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), val(sample_treatment_col), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.translated_mRNA.anota2seq.results.tsv"), emit: translated_mrna 
    tuple val(meta), path("*.total_mRNA.anota2seq.results.tsv")     , emit: total_mrna 
    tuple val(meta), path("*.translation.anota2seq.results.tsv")    , emit: translation
    tuple val(meta), path("*.buffering.anota2seq.results.tsv")      , emit: buffering
    tuple val(meta), path("*.mRNA_abundance.anota2seq.results.tsv") , emit: mrna_abundance
    tuple val(meta), path("*.Anota2seqDataSet.rds")                 , emit: rdata
    tuple val(meta), path("*.anota2seq.fold_change.png")            , emit: fold_change_plot
    tuple val(meta), path("*.R_sessionInfo.log")                    , emit: session_info
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'anota2seqrun.r'
}
