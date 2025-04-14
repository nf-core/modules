process DEEPBGC_DOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    output:
    path "deepbgc_db/", emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    export DEEPBGC_DOWNLOADS_DIR='./deepbgc_db'

    deepbgc \\
        download

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(echo \$(deepbgc info 2>&1 /dev/null/ | grep 'version' | cut -d " " -f3) )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p deepbgc_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepbgc: \$(echo \$(deepbgc info 2>&1 /dev/null/ | grep 'version' | cut -d " " -f3) )
    END_VERSIONS
    """
}
