process SAMSHEE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_samshee:6ef3bad21377d60e' :
        'community.wave.seqera.io/library/pip_samshee:9a11bec0d7c8df27' }"

    input:
    tuple val(meta), path(samplesheet)
    path(file_schema_validator)

    output:
    tuple val(meta), path("*_formatted.csv"), emit: samplesheet
    path "versions.yml"                     , emit: versions

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
