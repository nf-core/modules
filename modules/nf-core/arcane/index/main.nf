process ARCANE_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcane:1.0.0--pyh106432d_0':
        'quay.io/biocontainers/arcane:1.0.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(reference), path(genes)
    val n_objects

    output:
    tuple val(meta), path("*.hash"),                                        emit: index
    tuple val(meta), path("*.info"),                                        emit: info
    tuple val("${task.process}"), val('arcane'), eval("arcane --version"),  topic: versions, emit: versions_arcane

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cpus = task.cpus as int

    // The presented thread allocation behavior is taken from the tools source code
    // https://gitlab.com/rahmannlab/arcane/-/blob/main/arcane/arcane/arcane_index_colors.py?ref_type=heads
    def subtables = Math.max([(cpus / 2) as int - 1, cpus - 3, 19].min(), 1)
    if ((subtables % 2) == 0) {
        subtables += 1
    }

    def read_threads = Math.ceil(subtables / 10) as int

    def split_threads = 2 * read_threads

    if ((subtables + read_threads + split_threads) >= cpus) {
        read_threads = 1
        split_threads = 2
    }
    """
    arcane index \\
        --threads-split ${split_threads} \\
        --threads-read ${read_threads} \\
        --subtables ${subtables} \\
        --index ${prefix} \\
        --ref ${reference} \\
        --nobjects ${n_objects} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hash
    touch ${prefix}.info
    """
}
