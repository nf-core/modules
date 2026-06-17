process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5ec418b3aee9c5fd90b21d70e197009da826ff9a39945bb28181fbdc0b42d0b/data'
        : 'community.wave.seqera.io/library/bakta_diamond:cd8df03e82945085'}"

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
