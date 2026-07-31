process LDAK_KVIKSTEP1 {
    tag "${meta.id}:${meta2.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c94e5424f08cbd7e0856ab3c3a9992b4080944a5d0c497ce0abcceb413db4a3e/data'
        : 'community.wave.seqera.io/library/ldak6_r-base:452828f72b3c9129'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotype_file), val(is_binary)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)

    output:
    tuple val(meta2), path("${prefix}.step1.root"), path("${prefix}.step1.loco.details"), path("${prefix}.step1.loco.prs"), emit: predictions
    tuple val(meta2), path("${prefix}.step1.effects"), emit: effects, optional: true
    tuple val(meta2), path("${prefix}.step1.log"), emit: log
    tuple val("${task.process}"), val("ldak6"), eval("ldak6 --version 2>&1 | grep -oP '(?<=^Version )[0-9.]+'"), emit: versions_ldak6, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def bfile_prefix = bed.baseName
    prefix = task.ext.prefix ?: meta2.id
    def covar_arg = quant_covariates_file ? "--covar ${quant_covariates_file}" : ""
    def factors_arg = cat_covariates_file ? "--factors ${cat_covariates_file}" : ""
    def binary_arg = is_binary ? "--binary YES" : ""
    """
    ldak6 --kvik-step1 ${prefix} \\
        --bfile ${bfile_prefix} \\
        --pheno ${phenotype_file} \\
        ${covar_arg} \\
        ${factors_arg} \\
        ${binary_arg} \\
        --max-threads ${task.cpus} \\
        ${args} \\
        2>&1 | tee ${prefix}.step1.log
    """

    stub:
    prefix = task.ext.prefix ?: meta2.id
    """
    touch ${prefix}.step1.root
    touch ${prefix}.step1.loco.details
    touch ${prefix}.step1.loco.prs
    touch ${prefix}.step1.effects
    touch ${prefix}.step1.log
    """
}
