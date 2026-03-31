process DDS_CREATEPROJECT {

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
        path  "versions.yml", emit: versions

    script:
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
    """
    touch output.log

      cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dds: \$(dds --version)
    END_VERSIONS
    """
}
