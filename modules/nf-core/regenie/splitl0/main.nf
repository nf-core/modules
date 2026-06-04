process REGENIE_SPLITL0 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a05bf71ea09adc5ebf9f0c656c9b326c0f16ba8e4966914972e58313469a466/data'
        : 'community.wave.seqera.io/library/regenie:4.1.2--5d361f9fcb2f85cf'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    tuple val(meta2), path(pheno)
    tuple val(meta3), path(covar)
    val bsize
    val n_jobs

    output:
    tuple val(meta), path("*.master"), emit: master
    tuple val(meta), path("*_job*.snplist"), emit: snplists
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('regenie'), eval('regenie --version 2>&1 | sed -n "1{s/^v//;s/\\.gz$//;p}"'), topic: versions, emit: versions_regenie


    script:
    def args = task.ext.args ?: ''
    def input_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: input_prefix
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
        --out ${prefix} \\
        --split-l0 ${prefix},${n_jobs}
    """

    stub:
    def input_prefix = plink_genotype_file.baseName
    def prefix = task.ext.prefix ?: input_prefix
    def job_count = n_jobs as Integer
    def snplist_lines = (1..job_count).collect { job -> "touch ${prefix}_job${job}.snplist" }.join('\n')
    def master_lines = (1..job_count).collect { job -> "${prefix}_job${job} ${prefix}_job${job}.snplist" }.join('\\n')
    """
    printf 'job snplist\\n${master_lines}\\n' > ${prefix}.master
    ${snplist_lines}
    touch ${prefix}.log
    """
}
