process SAMSHEE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/pip_samshee:733e11f3377fc2e3' :
        'community.wave.seqera.io/library/pip_samshee:733e11f3377fc2e3' }"

    input:
    tuple val(meta), path(samplesheet)
    val(validator_schema)              //optional

    output:
    // Module is meant to stop the pipeline if validation fails
    tuple val(meta), path("*_formatted.csv"), emit: samplesheet
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def arg_validator_schema = validator_schema ? "--schema ${validator_schema}" : ""
    """
    # Run validation command and capture output
    python -m samshee $samplesheet \
    $arg_validator_schema \
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
