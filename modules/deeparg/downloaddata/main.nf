def VERSION='1.0.2'

process DEEPARG_DOWNLOADDATA {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::deeparg=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeparg:1.0.2--pyhdfd78af_1' :
        'quay.io/biocontainers/deeparg:1.0.2--pyhdfd78af_1' }"

    input:

    output:
    path "db/"          , emit: db
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    deeparg \\
        download_data \\
        $args \\
        -o db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeparg: $VERSION
    END_VERSIONS
    """
}
