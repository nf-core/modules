process HACKGAP_COUNT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hackgap:1.0.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/hackgap:1.0.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.hash"), emit: hash
    tuple val(meta), path("*.info"), emit: info
    tuple val("${task.process}"), val('hackgap'), eval("hackgap --version"), topic: versions, emit: versions_hackgap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cpus = task.cpus as int

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

    if (cpus < 4) {
        subtables = 1
        read_threads = 1
        split_threads = 1
    }
    """
    hackgap count \\
        ${args} \\
        --threads-split ${split_threads} \\
        --threads-read ${read_threads} \\
        --subtables ${subtables} \\
        --files ${reads} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hash
    touch ${prefix}.info
    """
}
