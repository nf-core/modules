process MEMOTE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/memote:0.17.0--pyhdfd78af_0' :
        'quay.io/biocontainers/memote:0.17.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(model)

    output:
    tuple val(meta), path("*.json.gz"), emit: json
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    memote run \\
        --filename ${prefix}.json.gz \\
        --ignore-git \\
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
    echo '{}' | gzip > ${prefix}.json.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        memote: \$(memote --version | sed 's/memote, version //')
    END_VERSIONS
    """
}
