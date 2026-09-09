process LEVIOSAM2_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leviosam2:0.4.2--h4ac6f70_0':
        'quay.io/biocontainers/leviosam2:0.4.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(fai)
    path(chain)

    output:
    tuple val(meta), path("*.clft"), emit: clft
    tuple val("${task.process}"), val('leviosam2'), eval("leviosam2 --version"), emit: versions_leviosam2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    leviosam2 \\
        index \\
        -c ${chain} \\
        -p ${prefix} \\
        -F ${fai}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clft
    """
}
