process GCTA_FASTGWA {
    tag "${meta.id}:${meta3.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(bed_pgen), path(bim_pvar), path(fam_psam)
    tuple val(meta2), path(sparse_grm_id), path(sparse_grm_sp)
    tuple val(meta3), path(phenotype_file)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)
    val mlm_exact

    output:
    tuple val(meta), val(meta3), path("${meta.id}_${meta3.id}.fastGWA"), emit: results
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def qcovar_arg = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_arg = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    def mpheno_arg = meta3.mpheno ? "--mpheno ${meta3.mpheno}" : ''
    def grm_arg = meta3.is_binary ? '' : "--grm-sparse ${meta2.id}"
    def genotype_suffix = bed_pgen.name.tokenize('.').last()
    def genotype_flag = genotype_suffix == 'pgen' ? '--pfile' : '--bfile'
    def genotype_prefix = meta.id
    def out = "${meta.id}_${meta3.id}"
    def extra_args = task.ext.args ?: ''
    def mode_arg = meta3.is_binary ? '--fastGWA-lr' : (mlm_exact ? '--fastGWA-mlm-exact' : '--fastGWA-mlm')

    """
    set -euo pipefail

    gcta \\
        ${genotype_flag} ${genotype_prefix} \\
        ${grm_arg} \\
        ${mode_arg} \\
        --pheno ${phenotype_file} \\
        ${qcovar_arg} \\
        ${covar_arg} \\
        ${mpheno_arg} \\
        --thread-num ${task.cpus} \\
        --out ${out} ${extra_args}
    """

    stub:
    """
    touch ${meta.id}_${meta3.id}.fastGWA
    """
}
