process LDAK_REML {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(phenotype_file), val(mpheno), val(prevalence)
    tuple val(meta2), path(grm_bin), path(grm_id), path(grm_details), path(grm_adjust)
    tuple val(meta3), path(keep_file)
    tuple val(meta4), path(quant_covariates_file)
    tuple val(meta5), path(cat_covariates_file)

    output:
    tuple val(meta), path("${prefix}.reml"), emit: reml_results
    tuple val(meta), path("${prefix}.coeff"), emit: coeff, optional: true
    tuple val(meta), path("${prefix}.combined"), emit: combined, optional: true
    tuple val(meta), path("${prefix}.cross"), emit: cross, optional: true
    tuple val(meta), path("${prefix}.indi.blp"), emit: indi_blp, optional: true
    tuple val(meta), path("${prefix}.indi.res"), emit: indi_res, optional: true
    tuple val(meta), path("${prefix}.progress"), emit: progress, optional: true
    tuple val(meta), path("${prefix}.share"), emit: share, optional: true
    tuple val(meta), path("${prefix}.vars"), emit: vars, optional: true
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def grm_prefix = grm_bin.name.replaceFirst(/\.grm\.bin$/, '')
    prefix = task.ext.prefix ?: "${meta.id}"
    def mpheno_arg = mpheno == null || (mpheno instanceof Collection && mpheno.isEmpty()) ? 1 : mpheno
    def keep_arg = keep_file ? "--keep ${keep_file}" : ''
    def quant_covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ''
    def cat_covar_arg = cat_covariates_file ? "--factors ${cat_covariates_file}" : ''
    def prevalence_arg = prevalence ? "--prevalence ${prevalence}" : ''

    """

    ldak6 --reml ${prefix} \\
        --pheno ${phenotype_file} \\
        --mpheno ${mpheno_arg} \\
        --grm ${grm_prefix} \\
        ${keep_arg} \\
        ${quant_covar_arg} \\
        ${cat_covar_arg} \\
        ${prevalence_arg} \\
        --max-threads ${task.cpus} \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reml
    touch ${prefix}.coeff
    touch ${prefix}.combined
    touch ${prefix}.cross
    touch ${prefix}.indi.blp
    touch ${prefix}.indi.res
    touch ${prefix}.progress
    touch ${prefix}.share
    touch ${prefix}.vars
    """
}
