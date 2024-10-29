process STIMULUS_ANALYSISDEFAULT {
    tag "$meta.id"
    label 'process_medium'

    // TODO: Add singularity support
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    conda "${moduleDir}/environment.yml"
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    path(data)
    tuple val(meta), path(experiment_config)
    path(model_config)
    path(weights)
    path(optimizer)
    path(metrics)
    path(model)

    output:
    tuple val(meta), path("performance_tune_train/"), path("performance_robustness/"), emit: analysis
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    stimulus-analysis-default \\
        $args \\
        -m ${model} \\
        -w ${weights} \\
        -me ${metrics} \\
        -ec ${experiment_config} \\
        -mc ${model_config} \\
        -d ${data} \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$( pip show stimulus-py | grep Version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def STIMULUS_VER = '0.0.9' // container not used in stub, change manually
    """
    mkdir performance_tune_train/
    mkdir performance_robustness/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
