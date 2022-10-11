process DESEQ2_DIFFERENTIAL {
    tag "$contrast_meta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-deseq2=1.32.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2%3A1.32.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-deseq2:1.32.0--r41h399db7b_0' }"

    input:
    path samplesheet
    path counts
    val contrast_meta

    output:
    tuple val(contrast_meta), path("*.deseq2.results.tsv")              , emit: results
    tuple val(contrast_meta), path("*.deseq2.dispersion.png")           , emit: dispersion_plot
    tuple val(contrast_meta), path("*.dds.rld.rds")                     , emit: rdata
    tuple val(contrast_meta), path("*.deseq2.sizefactors.tsv")          , emit: size_factors
    tuple val(contrast_meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    tuple val(contrast_meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    tuple val(contrast_meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    tuple val(contrast_meta), path("*.R_sessionInfo.log")               , emit: session_info
    tuple val(contrast_meta), path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_de.R'
}
