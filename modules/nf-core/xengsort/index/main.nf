process XENGSORT_INDEX {
    tag "$host_fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xengsort:2.0.5--pyhdfd78af_0':
        'quay.io/biocontainers/xengsort:2.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(host_fasta, stageAs: "host/*")
    tuple val(meta1), path(graft_fasta, stageAs: "graft/*")
    val nobjects
    val mask

    output:
    tuple val(meta1), path("*.hash")          , emit: hash
    tuple val(meta1), path("*.info")          , emit: info
    tuple val("${task.process}"), val('xengsort'), eval("xengsort --version"), topic: versions, emit: versions_xengsort

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    xengsort \\
        index \\
        $args \\
        --index ${prefix} \\
        --host ${host_fasta} \\
        --graft ${graft_fasta} \\
        --nobjects ${nobjects} \\
        --mask '${mask}'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.info
    touch ${prefix}.hash
    """
}
