process CLEANIFIER_INDEX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cleanifier:1.3.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/cleanifier:1.3.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reference)
    val n_objects

    output:
    tuple val(meta), path("*.{filter,hash}"), emit: index
    tuple val(meta), path("*.info"), emit: info
    tuple val("${task.process}"), val('cleanifier'), eval("cleanifier --version"), emit: versions_cleanifier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def threads_read = Math.max(1, (task.cpus / 4) as int)
    def threads_split = Math.max(1, (task.cpus - threads_read) as int)
    """
    cleanifier index \\
        --index ${prefix} \\
        --files ${reference} \\
        -n ${n_objects} \\
        --threads-read ${threads_read} \\
        --threads-split ${threads_split} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--exact") ? "hash" : "filter"
    """
    touch ${prefix}.${suffix}
    touch ${prefix}.info
    """
}
