process MINDAGAP_MINDAGAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mindagap:0.0.2--pyhdfd78af_1' :
    'biocontainers/mindagap:0.0.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(panorama)

    output:
    tuple val(meta), path("*.{tif,tiff}"), emit: tiff
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mindagap.py \\
        $panorama \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mindagap: \$(mindagap.py test -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${panorama.baseName}_gridfilled.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mindagap: \$(mindagap.py test -v)
    END_VERSIONS
    """
}
