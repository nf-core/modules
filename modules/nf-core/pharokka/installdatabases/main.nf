process PHAROKKA_INSTALLDATABASES {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.5.1--pyhdfd78af_0':
        'biocontainers/pharokka:1.5.1--pyhdfd78af_0' }"

    output:
    path("pharokka_db/")    , emit: pharokka_db
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    install_databases.py \\
        --outdir pharokka_db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p pharokka_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """
}
