process REGENIE_STEP1 {
    tag "${meta.id}:${meta2.pheno_col}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/regenie:4.1.2--7c121fb4ecd57890'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(pheno)
    tuple val(meta3), path(covar)
    val(bsize)

    output:
    tuple val(meta2), path("*_pred.list"), path("*.loco.gz"), emit: predictions
    tuple val(meta2), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval("regenie --version 2>&1 | head -n 1"), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def binary_arg = meta2.is_binary ? '--bt' : ''
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def pheno_col = meta2.pheno_col
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pred.list
    printf '' | gzip > ${prefix}_1.loco.gz
    touch ${prefix}.log
    """
}
