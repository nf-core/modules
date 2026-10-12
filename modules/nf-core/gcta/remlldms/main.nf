process GCTA_REMLLDMS {
    tag "${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(mgrm_file), path(grm_files)
    tuple val(meta2), path(phenotypes_file)
    tuple val(meta3), path(quant_covariates_file)
    tuple val(meta4), path(cat_covariates_file)

    output:
    tuple val(meta), path("*.hsq"), emit: reml_results
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def qcovar_param = quant_covariates_file ? "--qcovar ${quant_covariates_file}" : ''
    def covar_param = cat_covariates_file ? "--covar ${cat_covariates_file}" : ''
    """
    gcta \\
        --reml \\
        --mgrm ${mgrm_file} \\
        --pheno ${phenotypes_file} \\
        ${qcovar_param} \\
        ${covar_param} \\
        --out "${prefix}" \\
        --thread-num ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.hsq"
    """
}
