process DEEPBGC_DOWNLOAD {
    label 'process_single'

    conda "bioconda::deepbgc=0.1.30"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deepbgc:0.1.30--pyhb7b1952_1':
        'biocontainers/deepbgc:0.1.30--pyhb7b1952_1' }"

    output:
    path "deepbgc_db/"  , emit: db
    path "versions.yml" , emit: versions

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
}
