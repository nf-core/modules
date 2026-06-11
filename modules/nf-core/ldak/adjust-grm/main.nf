process LDAK_ADJUST_GRM {
    tag "${meta.id}:${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ldak6_r-base:e8012f1be139af77' :
        'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129' }"

    input:
    tuple val(meta), path(combined_grm_bin), path(combined_grm_id), path(combined_grm_details), path(combined_grm_adjust)
    tuple val(meta2), path(phenotype_file)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)
    val(mpheno)

    output:
    tuple val(meta), val(meta2), path("${meta.id}.grm.bin"), path("${meta.id}.grm.id"), path("${meta.id}.grm.details"), path("${meta.id}.grm.adjust"), path("${meta.id}.grm.root"), emit: adjusted_grm
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | head -n 1"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // --adjust-grm supports --covar but not --factors; categorical covariates are ignored here.
    covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ''
    grm_prefix = meta.grm_prefix
    mpheno_arg = (mpheno == null || (mpheno instanceof Collection && mpheno.isEmpty())) ? 1 : mpheno
    extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    ldak6 \
        --adjust-grm ${meta.id} \
        --grm ${grm_prefix} \
        --pheno ${phenotype_file} \
        --mpheno ${mpheno_arg} \
        ${covar_arg} \
        --max-threads ${task.cpus} ${extra_args}
    """

    stub:
    """
    touch ${meta.id}.grm.bin
    touch ${meta.id}.grm.id
    touch ${meta.id}.grm.details
    touch ${meta.id}.grm.adjust
    touch ${meta.id}.grm.root
    """
}
