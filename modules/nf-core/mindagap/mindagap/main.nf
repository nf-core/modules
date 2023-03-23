process MINDAGAP_MINDAGAP {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mindagap=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mindagap:0.0.2--pyhdfd78af_0' :
    'quay.io/biocontainers/mindagap:0.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.tiff"), emit: tiff
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mindagap.py \\
        $args \\
        $tiff \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mindagap: \$(mindagap.py test -v)
    END_VERSIONS
    """
}
