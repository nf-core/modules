process ASHLAR {
    tag '$meta.id'
    label 'process_single'

    conda "bioconda::ashlar=1.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ashlar:1.17.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/ashlar:1.17.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(file_in)

    output:
    tuple val(meta), path("*.ome.tif")      ,   emit: tif
    path "versions.yml"                     ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    ashlar \\
        $file_in \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(echo \$(ashlar --version 2>&1) | sed 's/^.*ashlar //; s/Using.*\$//' )
    END_VERSIONS
    """
}
