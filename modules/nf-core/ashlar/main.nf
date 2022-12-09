process ASHLAR {
    tag '$meta.id'
    label 'process_single'

    conda (params.enable_conda ? "bioconda::ashlar=1.17.0" : null)
    // TODO nf-core: need to have this added
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ashlar':
        'quay.io/biocontainers/ashlar' }"

    input:
    tuple val(meta), val(file_in)

    output:
    tuple val(meta), path("*.ome.tif")  ,   emit: ashlar_tif
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
        ashlar: \$(echo \$(ashlar --version 2>&1) | sed 's/^.*ashlar //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
