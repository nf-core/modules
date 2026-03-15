process REGENIE_STEP2 {
    tag "${meta.id}:${meta2.id}"
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

    output:
    tuple val(meta), val(meta2), path("*.regenie.gz"), emit: results
    tuple val(meta), val(meta2), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval("regenie --version 2>&1 | head -n 1"), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def binary_arg = meta2.is_binary ? '--bt' : ''
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def genotype_prefix = plink_genotype_file.baseName
    def base_prefix = meta.chr_prefix ?: meta.id
    def prefix = task.ext.prefix ?: "${base_prefix}.regenie.step2"
    def bsize = task.ext.bsize ?: 400
    def test_model = task.ext.regenie_test ?: 'additive'

    """
    regenie \\
        --step 2 \\
        ${genotype_flag} ${genotype_prefix} \\
        --phenoFile ${pheno} \\
        --phenoColList ${meta2.id} \\
        --pred ${predictions} \\
        ${covar_arg} \\
        ${binary_arg} \\
        --test ${test_model} \\
        --bsize ${bsize} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def base_prefix = meta.chr_prefix ?: meta.id
    def prefix = task.ext.prefix ?: "${base_prefix}.regenie.step2"
    """
    printf '' | gzip > ${prefix}_${meta2.id}.regenie.gz
    touch ${prefix}.log
    """
}
