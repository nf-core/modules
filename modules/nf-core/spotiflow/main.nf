process SPOTIFLOW {
    tag "$meta.id"
    label 'process_medium'
    container "docker.io/miguelib/spotiflow:0.0.3"

    input:
    tuple val(meta), path(image_2d)

    output:
    tuple val(meta), path("*.csv"), emit: spots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.0.3"

    """
    spotiflow-predict \\
        ${image_2d} \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spotiflow: $VERSION
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.0.3"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spotiflow: $VERSION
    END_VERSIONS
    """
}
