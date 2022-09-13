process EIDO_VALIDATE {
    tag '$samplesheet'
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::eido=0.1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/eido/0.1.9_cv1/eido_0.1.9_cv1.sif' :
        'biocontainers/eido:0.1.9_cv1' }"

    input:
    path samplesheet
    path schema

    output:
    path "versions.yml"           , emit: versions
    path "*.log"         , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "validation"
    """
    eido validate $args --st-index sample $samplesheet -s $schema -e > validation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' ))
    END_VERSIONS
    """
}
