process SAMSHEE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/pip_samshee:733e11f3377fc2e3' :
        'community.wave.seqera.io/library/pip_samshee:733e11f3377fc2e3' }"

    input:
    tuple val(meta), path(samplesheet)
    val(json_schema_validator)  // optional
    val(name_schema_validator)  // optional
    path(file_schema_validator) // optional

    output:
    tuple val(meta), path("*_formatted.csv"), emit: samplesheet
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def arg_json_schema_validator = json_schema_validator ? "--schema '${json_schema_validator}'"                         : ""
    def arg_name_schema_validator = name_schema_validator ? "--schema '${name_schema_validator}'"                         : ""
    def arg_file_schema_validator = file_schema_validator ? "--schema '{\"\$ref\": \"file:${file_schema_validator}\"}'" : ""
    def arg_v1_schema             = params.v1_schema      ? "--output-format sectioned"                                   : ""
    def args = task.ext.args ?: ""
    """
    # Run validation command and capture output
    python -m samshee $samplesheet \
    $arg_json_schema_validator \
    $arg_name_schema_validator \
    $arg_file_schema_validator \
    $arg_v1_schema \
    $args \
    > ${samplesheet.baseName}_formatted.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samshee: \$( python -m pip show --version samshee | grep "Version" | sed -e "s/Version: //g" )
        python: \$( python --version | sed -e "s/Python //g" )
    END_VERSIONS
    """

    stub:
    """
    touch ${samplesheet.baseName}_formatted.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samshee: \$( python -m pip show --version samshee | grep "Version" | sed -e "s/Version: //g" )
        python: \$( python --version | sed -e "s/Python //g" )
    END_VERSIONS
    """
}
