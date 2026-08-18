process REGENIE_RUNL0 {
    tag "${meta.id}_${job_number}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a05bf71ea09adc5ebf9f0c656c9b326c0f16ba8e4966914972e58313469a466/data'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(master), path(snplist), val(job_number)
    tuple val(meta3), path(pheno)
    tuple val(meta4), path(covar)
    val bsize

    output:
    tuple val(meta), path("*_l0_Y*"), emit: l0_predictions
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval('regenie --version 2>&1 | sed -n "1{s/^v//;s/\\.gz$//;p}"'), topic: versions, emit: versions_regenie

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: input_prefix
    def run_prefix = "${prefix}_job${job_number}"
    def genotype_flag = plink_genotype_file.name.endsWith('.pgen') ? '--pgen' : '--bed'
    def covar_arg = covar ? "--covarFile ${covar}" : ''
    def bsize_arg = bsize ?: 1000
    """
    regenie \\
        --step 1 \\
        ${genotype_flag} ${input_prefix} \\
        --phenoFile ${pheno} \\
        ${covar_arg} \\
        --bsize ${bsize_arg} \\
        --gz \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${run_prefix} \\
        --run-l0 ${master},${job_number}
    """

    stub:
    def input_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: input_prefix
    def run_prefix = "${prefix}_job${job_number}"
    """
    touch ${run_prefix}_l0_Y1
    touch ${run_prefix}.log
    """
}
