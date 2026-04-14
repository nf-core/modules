process REGENIE_STEP2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/regenie:4.1.2--7c121fb4ecd57890'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(pheno)
    tuple val(meta3), path(predictions), path(loco_files)
    tuple val(meta4), path(covar)
    val pheno_col
    val is_binary
    val bsize
    val test_model

    output:
    tuple val(meta), path("*.regenie.gz"), emit: results
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval("regenie --version 2>&1 | head -n 1"), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def binary_arg = is_binary ? '--bt' : ''
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def prefix = plink_genotype_file.baseName
    def bsize_arg = bsize ?: 400
    def test_model_arg = test_model ?: 'additive'

    """
    regenie \\
        --step 2 \\
        ${genotype_flag} ${prefix} \\
        --phenoFile ${pheno} \\
        --phenoColList ${pheno_col} \\
        --pred ${predictions} \\
        ${covar_arg} \\
        ${binary_arg} \\
        --test ${test_model_arg} \\
        --bsize ${bsize_arg} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = plink_genotype_file.baseName
    """
    printf '' | gzip > ${prefix}_${pheno_col}.regenie.gz
    touch ${prefix}.log
    """
}
