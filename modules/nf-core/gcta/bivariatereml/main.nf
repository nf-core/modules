process GCTA_BIVARIATEREML {
    tag "bivariate_reml_${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(phenotype_file)
    tuple val(meta2), path(grm_id), path(grm_bin), path(grm_n_bin)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)

    output:
    tuple val(meta), path("*.hsq"), emit: bivariate_results
    tuple val(meta), path("*.log"), emit: log_file
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | grep 'version v' | tr -s ' ' | cut -d' ' -f3 | sed 's/^v//'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pair_id = prefix
    def qcovar_param = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_param = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    def extra_args = task.ext.args ?: ''
    def expected_grm_basename = grm_id.name.replaceFirst(/\.grm\.id$/, '')
    if (meta2.id != expected_grm_basename) {
        throw new IllegalArgumentException("GCTA_BIVARIATEREML contract violation: meta2.id '${meta2.id}' must match GRM basename '${expected_grm_basename}'")
    }

    """
    set -euo pipefail

    gcta \\
        --reml-bivar 1 2 \\
        --grm ${meta2.id} \\
        --pheno "${phenotype_file}" \\
        ${qcovar_param} \\
        ${covar_param} \\
        --out "${pair_id}" \\
        --thread-num ${task.cpus} ${extra_args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pair_id = prefix
    """
    touch "${pair_id}.hsq"
    touch "${pair_id}.log"
    """
}
