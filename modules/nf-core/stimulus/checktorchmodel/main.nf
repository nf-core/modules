process STIMULUS_CHECKTORCHMODEL {
    tag "$experiment_config - $original_csv"
    label 'process_medium'

    // TODO freeze to Wave
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    path(original_csv)
    path(experiment_config)
    path(model)
    path(ray_tune_config)
    path(initial_weights)

    output:
    tuple val(meta), path("*_modelcheck.log"), emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: model.replaceFirst(/\.py/, "")
    """
    stimulus-check-model \
        -d ${original_csv} \
        -m ${model} \
        -e ${experiment_config} \
        -c ${ray_tune_config} \
        --initial_weights ${initial_weights} \
        --gpus ${task.accelerator.request} \
        --cpus ${task.cpus} \
        --memory "${task.memory}" \
        $args > ${prefix}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$(pip show stimulus-py | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def STIMULUS_VER = '0.0.9' // container not used in stub, change manually
    """
    touch ${meta.id}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
