process LDAK_PCGC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(phenotype_file)
    tuple val(meta2), path(grm_bin), path(grm_id), path(grm_details), path(grm_adjust), path(grm_root)
    tuple val(meta3), path(keep_file)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)
    val(mpheno)

    output:
    tuple val(meta), path("${meta.id}.pcgc"), emit: pcgc_results
    tuple val(meta), path("${meta.id}.pcgc.marginal"), emit: pcgc_marginal, optional: true
    tuple val(meta), path("${meta.id}.coeff"), emit: coeff, optional: true
    tuple val(meta), path("${meta.id}.combined"), emit: combined, optional: true
    tuple val(meta), path("${meta.id}.cross"), emit: cross, optional: true
    tuple val(meta), path("${meta.id}.progress"), emit: progress, optional: true
    tuple val(meta), path("${meta.id}.share"), emit: share, optional: true
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def keep_arg = keep_file ? "--keep ${keep_file}" : ''
    def quant_covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ''
    def cat_covar_arg = cat_covariates_file ? "--factors ${cat_covariates_file}" : ''
    def prevalence_arg = (meta.is_binary && meta.prevalence) ? "--prevalence ${meta.prevalence}" : ''
    def grm_prefix = meta2.grm_prefix
    def mpheno_arg = (mpheno == null || (mpheno instanceof Collection && mpheno.isEmpty())) ? 1 : mpheno
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    ldak6 --pcgc ${meta.id} \\
        --pheno ${phenotype_file} \\
        --mpheno ${mpheno_arg} \\
        --grm ${grm_prefix} \\
        ${keep_arg} \\
        ${quant_covar_arg} \\
        ${cat_covar_arg} \\
        ${prevalence_arg} \\
        --max-threads ${task.cpus} ${extra_args}
    """

    stub:
    """
    touch ${meta.id}.pcgc
    touch ${meta.id}.pcgc.marginal
    touch ${meta.id}.coeff
    touch ${meta.id}.combined
    touch ${meta.id}.cross
    touch ${meta.id}.progress
    touch ${meta.id}.share
    """
}
