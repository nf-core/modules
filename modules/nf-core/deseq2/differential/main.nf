process DESEQ2_DIFFERENTIAL {
    tag "$meta"
    label 'process_medium'

    conda "bioconda::bioconductor-deseq2=1.34.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.34.0--r41hc247a5b_3' :
        'biocontainers/bioconductor-deseq2:1.34.0--r41hc247a5b_3' }"

    input:
    tuple val(meta), val(contrast_variable), val(reference), val(target)
    tuple val(meta2), path(samplesheet), path(counts)
    tuple val(control_genes_meta), path(control_genes_file)

    output:
    tuple val(meta), path("*.deseq2.results.tsv")              , emit: results
    tuple val(meta), path("*.deseq2.dispersion.png")           , emit: dispersion_plot
    tuple val(meta), path("*.dds.rld.rds")                     , emit: rdata
    tuple val(meta), path("*.deseq2.sizefactors.tsv")          , emit: size_factors
    tuple val(meta), path("*.normalised_counts.tsv")           , emit: normalised_counts
    tuple val(meta), path("*.rlog.tsv")                        , optional: true, emit: rlog_counts
    tuple val(meta), path("*.vst.tsv")                         , optional: true, emit: vst_counts
    tuple val(meta), path("*.R_sessionInfo.log")               , emit: session_info
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_de.R'
}
