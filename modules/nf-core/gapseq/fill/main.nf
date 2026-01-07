process GAPSEQ_FILL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gapseq:1.4.0--hdfd78af_0' :
        'biocontainers/gapseq:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(draft), path(medium)

    output:
    tuple val(meta), path("*.RDS")  , emit: filled
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def medium_arg = medium ? "-m $medium" : ''
    """
    gapseq \\
        fill \\
        -p ${task.cpus} \\
        $medium_arg \\
        $args \\
        $draft

    # Rename output file
    mv *-filled.RDS ${prefix}-filled.RDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gapseq: \$(gapseq --version 2>&1 | sed 's/gapseq //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-filled.RDS
    """
}
