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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    coptr \\
        extract \\
        $args \\
        bamfolder/ \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "${task.process}":
            coptr: \$(coptr |& sed -E '11!d ; s/CoPTR.*?\\(v(.*?)\\).*/\\1/')
    END_VERSIONS
    """
}
