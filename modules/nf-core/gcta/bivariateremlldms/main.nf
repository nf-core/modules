process GCTA_BIVARIATEREMLLDMS {
    tag "bivariate_reml_ldms_${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(phenotype_file)
    tuple val(meta2), path(mgrm_file), path(grm_files)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)

    output:
    tuple val(meta), path("${meta.id}.hsq"), emit: bivariate_results
    tuple val(meta), path("${meta.id}.log"), emit: log_file
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pair_id = meta.id
    def qcovar_param = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_param = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    gcta \\
        --reml-bivar 1 2 \\
        --mgrm ${mgrm_file} \\
        --pheno "${phenotype_file}" \\
        ${qcovar_param} \\
        ${covar_param} \\
        --reml-bivar-no-constrain \\
        --reml-maxit 500 \\
        --out "${pair_id}" \\
        --thread-num ${task.cpus} ${extra_args}
    """

    stub:
    """
    touch "${meta.id}.hsq"
    touch "${meta.id}.log"
    """
}
