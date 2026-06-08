process MEMOTE_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/memote:0.17.0--pyhdfd78af_0' :
        'quay.io/biocontainers/memote:0.17.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(model)

    output:
    tuple val(meta), path("*.html"), emit: report , topic: report
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=\${PWD}
    export COBRA_SOLVER=glpk_exact

    memote report snapshot \\
        --filename ${prefix}.html \\
        ${args} \\
        ${model}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        memote: \$(memote --version | sed 's/memote, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        memote: 0.17.0
    END_VERSIONS
    """
}
