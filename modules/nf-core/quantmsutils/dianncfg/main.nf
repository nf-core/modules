process QUANTMSUTILS_DIANNCFG {
    tag "$meta.id"
    label 'process_tiny'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.23--pyh7e72e81_0' :
        'biocontainers/quantms-utils:0.0.23--pyh7e72e81_0' }"

    input:
    tuple val(meta), val(enzyme), val(fixed_modifications), val(variable_modifications)

    output:
    tuple val(meta), path("diann_config.cfg"), emit: diann_cfg
    path "*.log", emit: log
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    echo "--cut K*,R*,!*P --fixed-mod Carbamidomethyl,57.021464,C --var-mod Oxidation,15.994915,M" > diann_config.cfg
    touch GENERATE_DIANN_CFG.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
