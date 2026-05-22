process GCTA_FASTGWA {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(bed_pgen), path(bim_pvar), path(fam_psam)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)
    tuple val(meta5), path(sparse_grm_files)

    output:
    tuple val(meta), path("*.fastGWA"), emit: results
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out = prefix
    def genotype_suffix = bed_pgen.name.tokenize('.').last()
    def genotype_flag = genotype_suffix == 'pgen' ? '--pfile' : '--bfile'
    def genotype_prefix = bed_pgen.baseName
    def grm_arg = sparse_grm_files ? "--grm-sparse ${meta5.id}" : ''
    def qcovar_arg = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_arg = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    """
    gcta \\
        ${genotype_flag} ${genotype_prefix} \\
        ${grm_arg} \\
        --fastGWA-mlm \\
        --pheno ${phenotype_file} \\
        ${qcovar_arg} \\
        ${covar_arg} \\
        --thread-num ${task.cpus} \\
        --out ${out} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out = prefix
    """
    touch ${out}.fastGWA
    touch ${out}.log
    """
}
