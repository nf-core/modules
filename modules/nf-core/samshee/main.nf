process SAMSHEE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4acdaccdf86f5cedadbb91fba6e818c8cec006c874a19b84e60b8e2660b19f4c/data' :
        'community.wave.seqera.io/library/samshee_python:38088c103ef2751a' }"
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    // See https://github.com/lit-regensburg/samshee/issues/16 to see if this has been resolved.
    // python -m pip show samshee | sed -n "s/Version: //p" currently returns 0.0.0

    input:
    tuple val(meta), path(samplesheet)
    path(file_schema_validator)

    output:
    tuple val(meta), path("*_formatted.csv"), emit: samplesheet
    tuple val("${task.process}"), val('samshee'), val('0.2.13'), emit: versions_samshee, topic: versions
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
