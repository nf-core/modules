process REGENIE_STEP1 {
    tag "${meta.id}:${pheno_col}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a05bf71ea09adc5ebf9f0c656c9b326c0f16ba8e4966914972e58313469a466/data'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(pheno), val(pheno_col), val(is_binary)
    tuple val(meta3), path(covar)
    val bsize

    output:
    tuple val(meta), path("*_pred.list"), val(pheno_col), val(is_binary), emit: predictions
    tuple val(meta), path("*.loco.gz"), val(pheno_col), val(is_binary), emit: loco
    tuple val(meta), path("*.log"), val(pheno_col), val(is_binary), emit: log
    tuple val("${task.process}"), val('regenie'), eval('regenie --version 2>&1 | sed -n "1{s/^v//;s/\\.gz$//;p}"'), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def binary_arg = is_binary ? '--bt' : ''
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def prefix = plink_genotype_file.baseName
    def bsize_arg = bsize ?: 1000
    """
    regenie \\
        --step 1 \\
        ${genotype_flag} ${prefix} \\
        --phenoFile ${pheno} \\
        --phenoColList ${pheno_col} \\
        ${covar_arg} \\
        ${binary_arg} \\
        --bsize ${bsize_arg} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = plink_genotype_file.baseName
    """
    touch ${prefix}_pred.list
    echo "" | gzip > ${prefix}_1.loco.gz
    touch ${prefix}.log
    """
}
