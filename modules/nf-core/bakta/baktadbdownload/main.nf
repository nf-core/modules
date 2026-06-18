process BAKTA_BAKTADBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/50b75335f6394ae83fd05f364db27ee2eb75f4170e3525bb2aea47ad717a9e64/data'
        : 'community.wave.seqera.io/library/bakta_diamond:7830b94718da4f96'}"

    output:
    path "db*", emit: db
    tuple val("${task.process}"), val('bakta'), eval("bakta --version 2>&1 | sed 's/.*bakta //'"), emit: versions_bakta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export FONTCONFIG_PATH=\$PWD/.fontconfig

    bakta_db \\
        download \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export FONTCONFIG_PATH=\$PWD/.fontconfig

    echo "bakta_db \\
        download \\
        ${args}"

    mkdir -p db
    touch db/version.json
    touch db/bakta.db
    """
}
