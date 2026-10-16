process QUANTMSUTILS_DIANNCFG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7caaeee055202206dc635a6f3d154eb229162e924eda945fe76c88bda497aba/data'
:         'community.wave.seqera.io/library/quantms-utils_sdrf-pipelines_python_setuptools:a179872b9d0fa97d' }"

    input:
    tuple val(meta), val(enzyme), val(fixed_modifications), val(variable_modifications)

    output:
    tuple val(meta), path("diann_config.cfg"), emit: diann_cfg
    path "*.log", emit: log
    tuple val("${task.process}"), val('quantms-utils'), eval("quantmsutilsc --version 2>&1 | sed -n 's/quantmsutils //p'"), emit: versions_quantmsutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    quantmsutilsc dianncfg \\
        --enzyme "${enzyme}" \\
        --fix_mod "${fixed_modifications}" \\
        --var_mod "${variable_modifications}" \\
        ${args} 2>&1 | tee GENERATE_DIANN_CFG.log
    """

    stub:
    """
    echo "--cut K*,R*,!*P --fixed-mod Carbamidomethyl,57.021464,C --var-mod Oxidation,15.994915,M" > diann_config.cfg
    touch GENERATE_DIANN_CFG.log
    """
}
