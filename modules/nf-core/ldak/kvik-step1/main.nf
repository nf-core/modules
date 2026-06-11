process LDAK_KVIK_STEP1 {
    tag "${meta.id}:${meta2.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)
    val(mpheno)

    output:
    tuple val(meta2), path("${meta2.id}.step1.root"), path("${meta2.id}.step1.loco.details"), path("${meta2.id}.step1.loco.prs"), emit: predictions
    tuple val(meta2), path("${meta2.id}.step1.effects"), emit: effects, optional: true
    tuple val(meta2), path("${meta2.id}.step1.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = meta2.id
    def covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ""
    def factors_arg = cat_covariates_file ? "--factors ${cat_covariates_file}" : ""
    def binary_arg = meta2.is_binary ? "--binary YES" : ""
    def mpheno_value = (mpheno == null || (mpheno instanceof Collection && mpheno.isEmpty())) ? 1 : mpheno

    """
    set -euo pipefail

    ldak6 --kvik-step1 ${prefix} \\
        --bfile ${bed.baseName} \\
        --pheno ${phenotype_file} \\
        ${covar_arg} \\
        ${factors_arg} \\
        ${binary_arg} \\
        --mpheno ${mpheno_value} \\
        --max-threads ${task.cpus} \\
        ${args} 2>&1 | tee ${prefix}.step1.log
    """

    stub:
    def prefix = meta2.id
    """
    touch ${prefix}.step1.root
    touch ${prefix}.step1.loco.details
    touch ${prefix}.step1.loco.prs
    touch ${prefix}.step1.effects
    touch ${prefix}.step1.log
    """
}
