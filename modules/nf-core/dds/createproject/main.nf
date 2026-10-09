process DDS_CREATEPROJECT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28fa2a1c9a43350fb7f7533b7663a6202040a06b4dfaf89c9b21de6e1667594/data'
:         'community.wave.seqera.io/library/pip_python_dds-cli:da61a5f0f8c75872' }"

    input:
    tuple val(title), val(description), val(pi)
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
