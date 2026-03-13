process GCTA_REML {
    tag "gcta_reml_${meta.id}_${meta2.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(phenotypes_file)
    tuple val(meta2), path(grm_id), path(grm_bin), path(grm_n_bin)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)

    output:
    tuple val(meta), path("${meta.id}.hsq"), emit: reml_results
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def qcovar_param = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_param = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    gcta \\
        --reml \\
        --grm ${meta2.id} \\
        --pheno ${phenotypes_file} \\
        ${qcovar_param} \\
        ${covar_param} \\
        --out "${meta.id}" \\
        --thread-num ${task.cpus} ${extra_args}
    """

    stub:
    """
    touch "${meta.id}.hsq"
    """
}
