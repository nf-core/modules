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
    tuple val(meta), path("*.Runtime_log.txt"), emit: Runtime_log
    tuple val(meta), path("*.Fragmisincorporation_plot.pdf"), emit: Fragmisincorporation_plot
    tuple val(meta), path("*.Length_plot.pdf"), emit: Length_plot
    tuple val(meta), path("*.misincorporation.txt"), emit: misincorporation
    tuple val(meta), path("*.5pCtoT_freq.txt"), emit: 5pCtoT_freq
    tuple val(meta), path("*.3pGtoA_freq.txt"), emit: 3pGtoA_freq
    tuple val(meta), path("*.dnacomp.txt"), emit: dnacomp
    tuple val(meta), path("*.lgdistribution.txt"), emit: lgdistribution
    tuple val(meta), path("*.Stats_out_MCMC_hist.pdf"), emit: Stats_out_MCMC_hist
    tuple val(meta), path("*.Stats_out_MCMC_iter.csv"), emit: Stats_out_MCMC_iter
    tuple val(meta), path("*.Stats_out_MCMC_trace.pdf"), emit: Stats_out_MCMC_trace
    tuple val(meta), path("*.Stats_out_MCMC_iter_summ_stat.csv"), emit: Stats_out_MCMC_iter_summ_stat
    tuple val(meta), path("*.Stats_out_post_pred.pdf"), emit: Stats_out_post_pred
    tuple val(meta), path("*.Stats_out_MCMC_correct_prob.csv"), emit: Stats_out_MCMC_correct_prob
    tuple val(meta), path("*.dnacomp_genome.txt"), emit: dnacomp_genome
    tuple val(meta), path("*.rescaled.bam"), emit: rescaled
    tuple val(meta), path("*.fasta"), optional: true, emit: fasta
    tuple val(meta), path("*/"), optional: true, emit: path
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
        ${getSoftwareName(task.process)}: \$(mapDamage --version)
    END_VERSIONS
    """
}
