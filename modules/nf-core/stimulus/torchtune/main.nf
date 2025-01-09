process STIMULUS_TORCHTUNE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/mathysgrapotte/stimulus-py:latest"

    input:
    tuple val(meta),  path(data), path(experiment_config)
    tuple path(model_config), path(model)

    output:
    tuple  val(meta), path("*-model.safetensors"), emit: weights
    tuple  val(meta), path("*-optimizer.pt"), emit: optimizer
    tuple  val(meta), path("*-metrics.csv"),  emit: metrics

    // output the debug files if they are present, making this an optional channel
    path("ray_results/*/debug/best_model_*.txt")                                , emit: bestmodel       , optional: true
    path("ray_results/*/debug/worker_with_seed_*/model.safetensors")            , emit: safetensors     , optional: true
    path("ray_results/*/debug/worker_with_seed_*/seeds.txt")                    , emit: debug           , optional: true

    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stimulus-tuning \
        -c ${ray_tune_config} \
        -m ${model} \
        -d ${data_csv} \
        -e ${experiment_config} \
        -o ${prefix}-model.safetensors \
        -bo ${prefix}-optimizer.pt \
        -bm ${prefix}-metrics.csv \
        -bc ${prefix}-config.json \
        --initial_weights ${initial_weights} \
        --gpus ${task.accelerator.request} \
        --cpus ${task.cpus} \
        --memory "${task.memory}" \
        $args
    """

    stub:
    def prefix = task.ext.prefix
    """
    touch ${prefix}-model.safetensors ${prefix}-optimizer.pt ${prefix}-metrics.csv ${prefix}-config.json
    cat <<-END_VERSIONS > versions.yml
    """
}
