process ASHLAR {
    tag '$meta.id'
    label 'process_single'

    container "quay.io/biocontainers/ashlar:1.17.0--pyh5e36f6f_0"

    input:
    tuple val(meta), path(file_in)

    output:
    tuple val(meta), path("*.tif")      ,   emit: ashlar_tif
    path "versions.yml"                 ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args_conf = task.ext.args ?: ''
    def args_meta = meta.args ?: ''

    """
    ashlar \\
        $file_in \\
        $args_conf \\
        $args_meta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(echo \$(ashlar --version 2>&1) | sed 's/^.*ashlar //; s/Using.*\$//' )
    END_VERSIONS
    """
}
