process DDS_CREATEPROJECT {
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/971f09f84bdc0d023c1ac45905d21bfe4f1a36e4c822edf75325e55b68639a05/data':
        'community.wave.seqera.io/library/python_pip_dds-cli:c4d21ad6e2eb6702' }"

    input:
    val title
    val description
    val pi

    output:
    path 'output.log', emit: log
    tuple val("${task.process}"), val('dds'), eval("dds --version"), topic: versions, emit: versions_dds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    printf '%s' $DDS_CLI_TOKEN > .dds_cli_token
    chmod 600 .dds_cli_token
    DDS_CLI_ENV="dev-instance" dds --token-path .dds_cli_token project create  --title "$title" --description "$description" --principal-investigator "$pi" > output.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dds: \$(dds --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch output.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dds: \$(dds --version)
    END_VERSIONS
    """
}
