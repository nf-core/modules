process LDAK_KVIK_STEP2 {
    tag "${meta.id}:${meta2.id}:${meta.chr != null ? meta.chr : 'all'}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(step1_root), path(step1_loco_details), path(step1_loco_prs)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)
    tuple val(meta6), path(keep_file)
    val(mpheno)

    output:
    tuple val(meta), val(meta2), path("${meta2.id}.step2.assoc"), emit: results
    tuple val(meta), val(meta2), path("${meta2.id}.step2.summaries"), emit: summaries, optional: true
    tuple val(meta), val(meta2), path("${meta2.id}.step2.pvalues"), emit: pvalues, optional: true
    tuple val(meta), val(meta2), path("${meta2.id}.step2.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = meta2.id
    def bfile_prefix = bed.baseName
    def covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ""
    def factors_arg = cat_covariates_file ? "--factors ${cat_covariates_file}" : ""
    def keep_arg = keep_file ? "--keep ${keep_file}" : ""
    def chr_arg = meta.chr != null ? "--chr ${meta.chr}" : ""
    def mpheno_value = (mpheno == null || (mpheno instanceof Collection && mpheno.isEmpty())) ? 1 : mpheno

    """
    set -euo pipefail

    ldak6 --kvik-step2 ${prefix} \\
        --bfile ${bfile_prefix} \\
        --pheno ${phenotype_file} \\
        ${covar_arg} \\
        ${factors_arg} \\
        --mpheno ${mpheno_value} \\
        ${keep_arg} \\
        ${chr_arg} \\
        --by-chr NO \\
        --max-threads ${task.cpus} \\
        ${args} 2>&1 | tee ${prefix}.step2.log
    """

    stub:
    def prefix = meta2.id
    """
    touch ${prefix}.step2.assoc
    touch ${prefix}.step2.summaries
    touch ${prefix}.step2.log
    """
}
