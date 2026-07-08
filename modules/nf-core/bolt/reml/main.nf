process BOLT_REML {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bolt-lmm:2.5--h15e0e67_0'
        : 'quay.io/biocontainers/bolt-lmm:2.5--h15e0e67_0'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(phenotypes_file)
    tuple val(meta3), path(covariates_file)

    output:
    tuple val(meta), path("${meta.id}.bolt.reml.log"), emit: reml_log
    tuple val("${task.process}"), val("bolt"), eval("bolt --version 2>&1 | grep -o 'BOLT-LMM, v[0-9][0-9.]*' | head -n1 | sed 's/BOLT-LMM, v//'"), emit: versions_bolt, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bed_list = bed instanceof Collection ? bed : [bed]
    def bim_list = bim instanceof Collection ? bim : [bim]
    def bed_bim_flags = [bed_list, bim_list]
        .transpose()
        .collect { bed_file, bim_file -> "--bed \"${bed_file}\" --bim \"${bim_file}\"" }
        .join(' ')
    def covar_file_arg = covariates_file ? "--covarFile \"${covariates_file}\"" : ''

    """
    bolt \\
        --reml \\
        ${bed_bim_flags} \\
        --fam "${fam}" \\
        --phenoFile "${phenotypes_file}" \\
        ${covar_file_arg} \\
        --numThreads ${task.cpus} \\
        ${args} \\
        &> "${meta.id}.bolt.reml.log"
    """

    stub:
    """
    touch "${meta.id}.bolt.reml.log"
    """
}
