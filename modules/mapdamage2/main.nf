include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MAPDAMAGE2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::mapdamage2=2.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mapdamage2:2.2.1--pyr40_0"
    } else {
        container "quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0"
    }

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("results_*/Runtime_log.txt")                    ,emit: runtime_log
    tuple val(meta), path("results_*/Fragmisincorporation_plot.pdf")      ,emit: fragmisincorporation_plot
    tuple val(meta), path("results_*/Length_plot.pdf"), optional: true    ,emit: length_plot
    tuple val(meta), path("results_*/misincorporation.txt")               ,emit: misincorporation
    tuple val(meta), path("results_*/lgdistribution.txt")                 ,emit: lgdistribution
    tuple val(meta), path("results_*/dnacomp.txt")                        ,emit: dnacomp
    tuple val(meta), path("results_*/Stats_out_MCMC_hist.pdf")            ,emit: stats_out_mcmc_hist
    tuple val(meta), path("results_*/Stats_out_MCMC_iter.csv")            ,emit: stats_out_mcmc_iter
    tuple val(meta), path("results_*/Stats_out_MCMC_trace.pdf")           ,emit: stats_out_mcmc_trace
    tuple val(meta), path("results_*/Stats_out_MCMC_iter_summ_stat.csv")  ,emit: stats_out_mcmc_iter_summ_stat
    tuple val(meta), path("results_*/Stats_out_MCMC_post_pred.pdf")       ,emit: stats_out_mcmc_post_pred
    tuple val(meta), path("results_*/Stats_out_MCMC_correct_prob.csv")    ,emit: stats_out_mcmc_correct_prob
    tuple val(meta), path("results_*/dnacomp_genome.csv")                 ,emit: dnacomp_genome
    tuple val(meta), path("results_*/rescaled.bam"), optional: true       ,emit: rescaled
    tuple val(meta), path("results_*/5pCtoT_freq.txt")                    ,emit: pctot_freq
    tuple val(meta), path("results_*/3pGtoA_freq.txt")                    ,emit: pgtoa_freq
    tuple val(meta), path("results_*/*.fasta"), optional: true            ,emit: fasta
    tuple val(meta), path("results_*/*/"), optional: true                 ,emit: path
    
    path "versions.yml",emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mapDamage \\
       $options.args \\
       -i $bam \\
       -r $fasta  
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(mapDamage --version))
    END_VERSIONS
    """
}
