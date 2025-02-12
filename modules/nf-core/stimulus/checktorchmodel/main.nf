process STIMULUS_CHECKTORCHMODEL {
    tag "$data_config - $data"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker.io/mathysgrapotte/stimulus-py:0.2.4.dev"

    input:
    tuple val(meta), path(data), path(data_config)
    tuple val(meta), path(model), path(model_config), path(initial_weights)

    output:
    path "*_modelcheck.log", emit: log
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: model.baseName.replaceFirst(/\.py/, "")
    def weights_arg = initial_weights ? "--initial_weights ${initial_weights}" : ""
    """
    # initialize Ray
    ray start --head --port=6379 --temp-dir /tmp/ray

    # wait for it to start
    sleep 10

    # run the model check
    stimulus-check-model \
        -d ${data} \
        -m ${model} \
        -e ${data_config} \
        -c ${model_config} \
        ${weights_arg} \
        --ray_results_dirpath ${workDir} \
        $args > ${prefix}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: \$(pip show stimulus-py | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args         = task.ext.args ?: ''
    prefix           = task.ext.prefix ?: model.baseName.replaceFirst(/\.py/, "")
    def STIMULUS_VER = '0.2.2' // container not used in stub, change manually
    """
    touch ${prefix}_modelcheck.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version | cut -d ' ' -f 2)
        Stimulus-py: ${STIMULUS_VER}
    END_VERSIONS
    """
}
