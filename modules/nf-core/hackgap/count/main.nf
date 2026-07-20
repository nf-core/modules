process HACKGAP_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hackgap:1.0.1--pyhdfd78af_0':
        'quay.io/biocontainers/hackgap:1.0.1--pyhdfd78af_0' }"

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
    def split_threads = Math.max(1, (task.cpus / 4) as int)
    def read_threads = Math.max(1, (task.cpus / 4) as int)
    def subtables_threads = Math.max(1, task.cpus - split_threads - read_threads)
    """
    hackgap count \\
        $args \\
        --threads-split ${split_threads} \\
        --threads-read ${read_threads} \\
        --subtables ${subtables_threads} \\
        --files ${reads} \\
        --out ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hash
    touch ${prefix}.info
    """
}
