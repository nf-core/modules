process DDS_CREATEPROJECT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/7352a40103ce629fbc72986922a68e8f9cc35630b89ee836a0c03dfd676f2341/data':
        'community.wave.seqera.io/library/pip_python_dds-cli:4a7daeb54b4c7a8c' }"

    input:
    val title
    val description
    val pi
    path token_file

    output:
    path 'output.log', emit: log
    tuple val("${task.process}"), val('dds'), eval("dds --version | sed 's/.*version //'"), topic: versions, emit: versions_dds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # dds requires 600 permissions on the token file; copy first since Nextflow stages inputs as symlinks
    cp $token_file token.conf
    chmod 600 token.conf
    dds --token-path token.conf project create --title "$title" --description "$description" --principal-investigator "$pi" > output.log
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch output.log
    """
}
