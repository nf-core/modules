process GCTA_REMLLDMS {
    tag "gcta_reml_ldms_${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(phenotypes_file)
    tuple val(meta2), path(mgrm_file), path(grm_files)
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
        --reml-no-constrain \\
        --mgrm ${mgrm_file} \\
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
