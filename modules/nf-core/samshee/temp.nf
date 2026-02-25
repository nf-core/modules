process SAMSHEE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/874139c488abfe0d65786bdbafea96a0768554c1a5e8e70a332708954b96dee7/data' :
        'community.wave.seqera.io/library/samshee_python:df66f39f919076ff' }"

    input:
    tuple val(meta), path(samplesheet)
    path(file_schema_validator)

    output:
    tuple val(meta), path("*_formatted.csv"), emit: samplesheet
    tuple val("${task.process}"), val('samshee'), eval('python -m pip show samshee | sed -n "s/Version: //p"'), emit: versions_samshee, topic: versions
    tuple val("${task.process}"), val('python'), eval('python --version | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def arg_file_schema_validator = file_schema_validator ? "--schema '{\"\$ref\": \"file:${file_schema_validator}\"}'"   : ""
    def args = task.ext.args ?: ""
    """
    # Run validation command and capture output
    python -m samshee $samplesheet \
    $args \
    $arg_file_schema_validator \
    > ${samplesheet.baseName}_formatted.csv
    """

    stub:
    """
    touch ${samplesheet.baseName}_formatted.csv
    """
}
