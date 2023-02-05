process MINDAGAP {
    tag "$meta.id"
    label 'process_low'

    container "rguerr/mindagap:latest"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.tiff"), emit: tiff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /mindagap/mindagap.py \\
        $tiff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mindagap: \$(python ./mindagap.py test -v))
    END_VERSIONS
    """
}
