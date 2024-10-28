process STIMULUS_CHECKTORCHMODEL {
    tag "$experiment_config - $original_csv"
    label 'process_medium'

    // TODO freeze to Wave
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    tuple val(meta), path(original_csv), path(experiment_config)
    tuple val(meta2), path(model), path(ray_tune_config), path(initial_weights)

    output:
    tuple val(meta2), path("${meta.id}_${meta2.id}_modelcheck.log"), emit: log
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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
        $args > ${meta.id}_${meta2.id}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$(pip show stimulus-py | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${meta.id}_${meta2.id}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$(pip show stimulus-py | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
