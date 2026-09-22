process MAPDAMAGE2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapdamage2:2.2.1--pyr40_0' :
        'quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("results_*/Runtime_log.txt")                  , emit: runtime_log
    tuple val(meta), path("results_*/Fragmisincorporation_plot.pdf")    , emit: fragmisincorporation_plot    , optional: true
    tuple val(meta), path("results_*/Length_plot.pdf")                  , emit: length_plot                  , optional: true
    tuple val(meta), path("results_*/misincorporation.txt")             , emit: misincorporation             , optional: true
    tuple val(meta), path("results_*/lgdistribution.txt")               , emit: lgdistribution               , optional: true
    tuple val(meta), path("results_*/dnacomp.txt")                      , emit: dnacomp                      , optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_hist.pdf")          , emit: stats_out_mcmc_hist          , optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_iter.csv")          , emit: stats_out_mcmc_iter          , optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_trace.pdf")         , emit: stats_out_mcmc_trace         , optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_iter_summ_stat.csv"), emit: stats_out_mcmc_iter_summ_stat, optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_post_pred.pdf")     , emit: stats_out_mcmc_post_pred     , optional: true
    tuple val(meta), path("results_*/Stats_out_MCMC_correct_prob.csv")  , emit: stats_out_mcmc_correct_prob  , optional: true
    tuple val(meta), path("results_*/dnacomp_genome.csv")               , emit: dnacomp_genome               , optional: true
    tuple val(meta), path("results_*/*rescaled.bam")                    , emit: rescaled                     , optional: true
    tuple val(meta), path("results_*/5pCtoT_freq.txt")                  , emit: pctot_freq                   , optional: true
    tuple val(meta), path("results_*/3pGtoA_freq.txt")                  , emit: pgtoa_freq                   , optional: true
    tuple val(meta), path("results_*/*.fasta")                          , emit: fasta                        , optional: true
    tuple val(meta), path("*/")                                         , emit: folder                       , optional: true
    tuple val("${task.process}"), val("mapdamage2"), eval("mapDamage --version"), topic: versions, emit: versions_mapdamage2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mapDamage \\
        ${args} \\
        -i ${bam} \\
        -r ${fasta}
    """

    stub:
    """
    mkdir -p results_${bam.baseName}
    cd results_${bam.baseName}
    touch Runtime_log.txt

    touch Fragmisincorporation_plot.pdf
    touch Length_plot.pdf

    touch misincorporation.txt
    touch lgdistribution.txt

    touch Stats_out_MCMC_{hist.pdf,iter.csv,trace.pdf,iter_summ_stat.csv,post_pred.pdf,correct_prob.csv}

    touch dnacomp_genome.csv
    touch dnacomp.txt

    touch 5pCtoT_freq.txt
    touch 3pGtoA_freq.txt
    """
}
