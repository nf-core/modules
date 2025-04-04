process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"


    container "ghcr.io/bcb-unl/run_dbcan_new:5.0.2"

    output:

    path "dbcan_db", emit: dbcan_db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '5.0.2'
    """

    run_dbcan database \\
        --db_dir dbcan_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '5.0.2'

    """
    mkdir  -p dbcan_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: $VERSION
    END_VERSIONS
    """

}
