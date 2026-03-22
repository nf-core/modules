process BOLT_REML {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bolt-lmm:2.5--h15e0e67_0' :
        'biocontainers/bolt-lmm:2.5--h15e0e67_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotypes_file)
    tuple val(meta3), path(covariates_file)

    output:
    tuple val(meta), path("${meta.id}.bolt.reml.log"), emit: reml_log
    tuple val("${task.process}"), val("bolt"), eval("bolt --version 2>&1 | head -n 1 | sed 's/^BOLT-LMM v//'"), emit: versions_bolt, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bed_list = bed instanceof Collection ? bed : [bed]
    def bim_list = bim instanceof Collection ? bim : [bim]
    def bed_bim_flags = [bed_list, bim_list]
        .transpose()
        .collect { bed_file, bim_file -> "--bed \"${bed_file}\" --bim \"${bim_file}\"" }
        .join(' ')
    def pheno_col_arg = meta2.pheno_col ? "--phenoCol=\"${meta2.pheno_col}\"" : ''
    def covariates_columns = meta3.covariates_columns ?: ''
    def covariates_cat_columns = meta3.covariates_cat_columns ?: ''
    def covariate_names = covariates_columns ? covariates_columns.split(',')*.trim().findAll { it } : []
    def categorical_covariate_names = covariates_cat_columns ? covariates_cat_columns.split(',')*.trim().findAll { it } : []
    def quantitative_covariate_names = covariate_names - categorical_covariate_names
    def covar_file_arg = (covariates_file && covariate_names) ? "--covarFile \"${covariates_file}\"" : ''
    def covar_cols_arg = categorical_covariate_names.collect { column_name -> "--covarCol=\"${column_name.replace('"', '\\"')}\"" }.join(' ')
    def qcovar_cols_arg = quantitative_covariate_names.collect { column_name -> "--qCovarCol=\"${column_name.replace('"', '\\"')}\"" }.join(' ')
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    bolt \\
        --reml \\
        ${bed_bim_flags} \\
        --fam "${fam}" \\
        --phenoFile "${phenotypes_file}" \\
        ${pheno_col_arg} \\
        ${covar_file_arg} \\
        ${covar_cols_arg} \\
        ${qcovar_cols_arg} \\
        --numThreads ${task.cpus} ${extra_args} \\
        &> "${meta.id}.bolt.reml.log"
    """

    stub:
    """
    touch "${meta.id}.bolt.reml.log"
    """
}
