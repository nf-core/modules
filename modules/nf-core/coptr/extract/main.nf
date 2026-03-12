process COPTR_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coptr:1.1.4--pyhdfd78af_3':
        'biocontainers/coptr:1.1.4--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(bam, stageAs: "bamfolder/*")

    output:
    tuple val(meta), path("*.pkl"), emit: coverage
    tuple val("${task.process}"), val('coptr'), eval("coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/'"), emit: versions_coptr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''

    """
    coptr \\
        extract \\
        $args \\
        bamfolder/ \\
        .

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pkl
    """
}
