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
    tuple val(meta), path("*.translated_mRNA.anota2seq.results.tsv")            , emit: translated_mrna, optional: true
    tuple val(meta), path("*.total_mRNA.anota2seq.results.tsv")                 , emit: total_mrna, optional: true
    tuple val(meta), path("*.translation.anota2seq.results.tsv")                , emit: translation, optional: true
    tuple val(meta), path("*.buffering.anota2seq.results.tsv")                  , emit: buffering, optional: true
    tuple val(meta), path("*.mRNA_abundance.anota2seq.results.tsv")             , emit: mrna_abundance, optional: true
    tuple val(meta), path("*.Anota2seqDataSet.rds")                             , emit: rdata, optional: true
    tuple val(meta), path("*.fold_change.png")                                  , emit: fold_change_plot, optional: true
    tuple val(meta), path("*.interaction_p_distribution.pdf")                   , emit: interaction_p_distribution_plot, optional: true
    tuple val(meta), path("*.residual_distribution_summary.jpeg")               , emit: residual_distribution_summary_plot, optional: true
    tuple val(meta), path("*.residual_vs_fitted.jpeg")                          , emit: residual_vs_fitted_plot, optional: true
    tuple val(meta), path("*.rvm_fit_for_all_contrasts_group.jpg")              , emit: rvm_fit_for_all_contrasts_group_plot, optional: true
    tuple val(meta), path("*.rvm_fit_for_interactions.jpg")                     , emit: rvm_fit_for_interactions_plot, optional: true
    tuple val(meta), path("*.rvm_fit_for_omnibus_group.jpg")                    , emit: rvm_fit_for_omnibus_group_plot, optional: true
    tuple val(meta), path("*.simulated_vs_obt_dfbetas_without_interaction.pdf") , emit: simulated_vs_obt_dfbetas_without_interaction_plot, optional: true
    tuple val(meta), path("*.R_sessionInfo.log")                                , emit: session_info
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/modules/nf-core/anota2seq/anota2seqrun/templates/anota2seqrun.r --output_prefix ${task.ext.prefix ?: meta.id} --sample_treatment_col ${sample_treatment_col} --reference_level ${reference} --target_level ${target} --sample_file ${samplesheet} --count_file ${counts}
    """
}
