process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0'
        : 'biocontainers/bakta:1.11.4--pyhdfd78af_0'}"

    output:
    path "db*", emit: db
    tuple val("${task.process}"), val('bakta'), eval("bakta --version 2>&1 | sed 's/bakta //'"), emit: versions_bakta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    bakta_db \\
        download \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    echo "bakta_db \\
        download \\
        ${args}"

    mkdir -p db
    touch db/version.json
    touch db/bakta.db
    """
}
