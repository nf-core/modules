params.deseq_gene_id_col = "gene_id"
params.deseq_sample_id_col = "experiment_accession"
params.deseq_test = "Wald"
params.deseq_fit_type = "parametric"
params.deseq_sf_type = 'ratio'
params.deseq_min_replicates_for_replace = 7
params.deseq_use_t = 'FALSE'
params.deseq_lfc_threshold = 0
params.deseq_alt_hypothesis = 'greaterAbs'
params.deseq_independent_filtering = 'TRUE'
params.deseq_p_adjust_method = 'BH'
params.deseq_alpha = 0.1
params.deseq_minmu = 0.5
params.deseq_write_normalised = 'TRUE'
params.deseq_write_variance_stabilised = 'TRUE'
params.deseq_vs_method = 'vst'

process DESEQ2_DIFFERENTIAL {
    tag "$meta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-deseq2=1.28.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-deseq2%3A1.28.0--r40h5f743cb_0' :
        'quay.io/biocontainers/bioconductor-deseq2:1.28.0--r40h5f743cb_0' }"

    input:
    path samplesheet
    path counts
    val meta

    output:
    tuple val(meta), path("deseq2.results.tsv")               , emit: results
    tuple val(meta), path("deseq2.dispersion.png")            , emit: dispersion_plot
    tuple val(meta), path("dds.rld.RData")                    , emit: rdata
    tuple val(meta), path("deseq2.sizefactors.tsv")           , emit: size_factors
    tuple val(meta), path("normalised_counts.tsv")            , emit: normalised_counts
    tuple val(meta), path("variance_stabilised_counts.tsv")   , emit: variance_stabilised_counts
    tuple val(meta), path("R_sessionInfo.log")                , emit: session_info

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'deseq_de.r'
}
