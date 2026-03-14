process REGENIE_STEP1 {
    tag "${meta.id}:${meta2.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/regenie:4.1.2--7c121fb4ecd57890'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(pheno)
    tuple val(meta3), path(covar)

    output:
    tuple val(meta2), path("*_pred.list"), path("*.loco.gz"), emit: predictions
    tuple val(meta2), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval("regenie --version 2>&1 | head -n 1"), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def binary_arg = meta2.is_binary ? '--bt' : ''
    def covar_file = covar instanceof List ? covar.find() : covar
    def covar_arg = covar_file ? "--covarFile ${covar_file}" : ''
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def genotype_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: "${meta2.id}.regenie.step1"
    def bsize = task.ext.bsize ?: 1000

    """
    regenie \\
        --step 1 \\
        ${genotype_flag} ${genotype_prefix} \\
        --phenoFile ${pheno} \\
        --phenoColList ${meta2.id} \\
        ${covar_arg} \\
        ${binary_arg} \\
        --bsize ${bsize} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta2.id}.regenie.step1"
    """
    touch ${prefix}_pred.list
    printf '' | gzip > ${prefix}_1.loco.gz
    touch ${prefix}.log
    """
}
