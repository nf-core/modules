process ARGNORM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/argnorm:0.5.0--pyhdfd78af_0':
        'biocontainers/argnorm:0.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_tsv)
    val(db)
    val(hamronized)
    val(tool)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.5.0'
    """
    argnorm \\
        $tool \\
        -i $input_tsv \\
        -o $prefix \\
        --db $db \\
        $hamronized

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        argnorm: $VERSION
    END_VERSIONS
    """
}
