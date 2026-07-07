process LEAFCUTTER_DIFFERENTIALSPLICING {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leafcutter:2.0.3--pyhd8ed1ab_0':
        'quay.io/biocontainers/leafcutter:2.0.3--pyhd8ed1ab_0' }"

    input:
    tuple val(meta), path(counts), path(groups)

    output:
    tuple val(meta), path("*_cluster_significance.txt"), emit: cluster_significance
    tuple val(meta), path("*_effect_sizes.txt")        , emit: effect_sizes
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('leafcutter'), val("2.0.3"), topic: versions, emit: versions_leafcutter

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export USER=\${USER:-nobody}

    leafcutter-ds \\
        $counts \\
        $groups \\
        --output_prefix ${prefix}_results \\
        --num_threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results_cluster_significance.txt
    touch ${prefix}_results_effect_sizes.txt
    """
}
